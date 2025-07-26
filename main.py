import tkinter as tk
from tkinter import filedialog
import pulp
from PIL import Image, ImageDraw


def parse_formula(formula):
    stack = [{}]
    i = 0
    n = len(formula)
    while i < n:
        if formula[i].isupper():
            elem = formula[i]
            i += 1
            if i < n and formula[i].islower():
                elem += formula[i]
                i += 1
            num = 0
            while i < n and formula[i].isdigit():
                num = num * 10 + int(formula[i])
                i += 1
            num = num if num > 0 else 1
            stack[-1][elem] = stack[-1].get(elem, 0) + num
        elif formula[i] in '([{':
            stack.append({})
            i += 1
        elif formula[i] in ')]}':
            i += 1
            num = 0
            while i < n and formula[i].isdigit():
                num = num * 10 + int(formula[i])
                i += 1
            num = num if num > 0 else 1
            current = stack.pop()
            for e in current:
                current[e] *= num
            for e in current:
                stack[-1][e] = stack[-1].get(e, 0) + current[e]
        else:
            i += 1  # skip invalid
    if len(stack) != 1:
        raise ValueError("Unbalanced parentheses in formula")
    return stack[0]


class StoichApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Chemical Stoichiometry Balancer")

        # Add notes about formula entry
        notes_text = """
Formula Entry Rules:
- Enter chemical formulas in linear notation, e.g., H2O, C6H12O6, (NH4)2SO4.
- Elements: Start with uppercase letter, optional lowercase (e.g., Na, Cl).
- Counts: Numbers after elements or groups (default 1 if omitted).
- Groups: Use parentheses (), brackets [], or braces {} for subgroups, followed by optional count.
- Nested groups are supported.
- No spaces; invalid characters are skipped, but avoid them to prevent errors.

Special Cases:
- Ensure parentheses are balanced; unbalanced will raise an error.
- For hydrates, use dot notation if needed, but parser treats '.' as invalid (skip), so use groups instead, e.g., CuSO4(H2O)5.
- Examples: Al2(SO4)3, [Co(NH3)6]Cl3, C6H5OH.
"""
        tk.Label(root, text="Notes on Formula Syntax:", font=("Arial", 12, "bold")).pack()
        notes = tk.Text(root, height=15, width=60, wrap=tk.WORD)
        notes.insert(tk.END, notes_text.strip())
        notes.config(state=tk.DISABLED)
        notes.pack()

        tk.Label(root, text="Number of reactants:").pack()
        self.num_entry = tk.Entry(root)
        self.num_entry.pack()
        tk.Button(root, text="Set Reactants", command=self.set_reactants).pack()

        self.reactant_entries = []
        self.target_entry = None
        self.result_text = None
        self.react_formulas = []
        self.target_form = ""
        self.moles = []
        self.leftovers = {}
        self.a_vars = []

    def set_reactants(self):
        try:
            num = int(self.num_entry.get())
        except ValueError:
            return
        for widget in self.reactant_entries:
            widget.destroy()
        self.reactant_entries = []
        for i in range(num):
            tk.Label(self.root, text=f"Reactant {i + 1} formula:").pack()
            e = tk.Entry(self.root)
            e.pack()
            self.reactant_entries.append(e)
        tk.Label(self.root, text="Target formula:").pack()
        self.target_entry = tk.Entry(self.root)
        self.target_entry.pack()
        tk.Button(self.root, text="Calculate", command=self.calculate).pack()
        if self.result_text:
            self.result_text.destroy()
        self.result_text = tk.Text(self.root, height=10, width=50)
        self.result_text.pack()
        tk.Button(self.root, text="Save as JPG", command=self.save_jpg).pack()

    def calculate(self):
        self.react_formulas = [e.get().strip() for e in self.reactant_entries]
        self.target_form = self.target_entry.get().strip()
        try:
            target = parse_formula(self.target_form)
            reactants = [parse_formula(f) for f in self.react_formulas]
        except Exception as e:
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, f"Error parsing formulas: {e}")
            return

        all_atoms = set(target.keys())
        for r in reactants:
            all_atoms.update(r.keys())

        target_atoms = list(target.keys())
        num_react = len(reactants)
        if num_react == 0:
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, "No reactants provided")
            return

        prob = pulp.LpProblem("Stoichiometry", pulp.LpMinimize)
        self.a_vars = [pulp.LpVariable(f"a_{i}", lowBound=0) for i in range(num_react)]
        num_atoms_ri = [sum(r.values()) for r in reactants]
        prob += pulp.lpSum(self.a_vars[i] * num_atoms_ri[i] for i in range(num_react))

        for atom in target_atoms:
            prob += pulp.lpSum(self.a_vars[i] * reactants[i].get(atom, 0) for i in range(num_react)) >= target[
                atom], f"Constraint_{atom}"

        status = prob.solve(pulp.PULP_CBC_CMD(msg=0))
        if status != pulp.LpStatusOptimal:
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, "No feasible solution found. Cannot form the target with given reactants.")
            return

        self.moles = [pulp.value(var) for var in self.a_vars]

        provided = {atom: 0.0 for atom in all_atoms}
        for i, r in enumerate(reactants):
            for atom, count in r.items():
                provided[atom] += self.moles[i] * count

        self.leftovers = {}
        for atom, prov in provided.items():
            req = target.get(atom, 0)
            excess = prov - req
            if excess > 1e-6:
                self.leftovers[atom] = excess

        self.result_text.delete(1.0, tk.END)
        text = "Theoretical moles needed for 1 mole of target:\n"
        for i, mol in enumerate(self.moles):
            text += f"{self.react_formulas[i]}: {mol:.4f}\n"
        text += "\nLeftOver atoms:\n"
        if not self.leftovers:
            text += "None\n"
        else:
            for atom, cnt in self.leftovers.items():
                text += f"{atom}: {cnt:.4f}\n"
        self.result_text.insert(tk.END, text)

    def save_jpg(self):
        if not self.moles:
            return
        filename = filedialog.asksaveasfilename(defaultextension=".jpg", filetypes=[("JPEG files", "*.jpg")])
        if not filename:
            return
        text = " + ".join(
            f"{self.moles[i]:.2f} {self.react_formulas[i]}" for i in range(len(self.moles)) if self.moles[i] > 1e-6)
        text += f" -> 1 {self.target_form}"
        if self.leftovers:
            text += " + leftovers: " + " + ".join(f"{cnt:.2f} {atom}" for atom, cnt in self.leftovers.items())
        else:
            text += " (no leftovers)"

        # Create image with rough size estimation
        width = len(text) * 8 + 20  # approximate
        img = Image.new('RGB', (width, 50), color='white')
        d = ImageDraw.Draw(img)
        d.text((10, 10), text, fill='black')
        img.save(filename)


if __name__ == "__main__":
    root = tk.Tk()
    app = StoichApp(root)
    root.mainloop()
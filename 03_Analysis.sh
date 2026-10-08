# 1. Get free energy (dG) 
openfe gather results/ -o absolute_free_energies.tsv

# 2. Get relative free energy (ddG) 
openfe gather results/ --report ddg -o relative_free_energies.tsv

# 3. Raw values of energy
openfe gather results/ --report raw -o raw_results.tsv

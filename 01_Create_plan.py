#!/usr/bin/env python
# coding: utf-8

# In[101]:


"""
Doncho's
This code is intended to prepare tranformations for RBFE 
# Requirements
1. A PDB structure of the protein
2. A SDF file with all the transformation
"""
# 0.1 Paths to files and variables

aligands = "tyk2_ligands.sdf"
pdb_prot = "tyk2_protein.pdb"
cpu_num = 14

# 0.2 Import openfe
import openfe


# In[102]:


# 1.0 Load chemical structures
from rdkit import Chem
supp = Chem.SDMolSupplier(aligands, removeHs=False)
ligands = [openfe.SmallMoleculeComponent.from_rdkit(mol) for mol in supp]

# 1.1 Add charges to molecules default uses 8 cpus
from openfe.protocols.openmm_utils.omm_settings import OpenFFPartialChargeSettings
from openfe.protocols.openmm_utils.charge_generation import bulk_assign_partial_charges
charge_settings = OpenFFPartialChargeSettings(partial_charge_method="am1bcc", off_toolkit_backend="ambertools")
charged_ligands = bulk_assign_partial_charges(
    molecules=ligands,
    overwrite=False,  
    method=charge_settings.partial_charge_method,
    toolkit_backend=charge_settings.off_toolkit_backend,
    generate_n_conformers=charge_settings.number_of_conformers,
    nagl_model=charge_settings.nagl_model,
    processors=cpu_num
)

# 2.1 Define the mapper for perturbations
mapper = openfe.KartografAtomMapper()
scorer = openfe.lomap_scorers.default_lomap_score
network_planner = openfe.ligand_network_planning.generate_minimal_spanning_network

#


# In[103]:


# 2.2 Create the network for perturbations
from pathlib import Path

ligand_network = network_planner(
    ligands=charged_ligands,
    mappers=[mapper],
    scorer=scorer
)

from openfe.utils.atommapping_network_plotting import plot_atommapping_network

plot_atommapping_network(ligand_network)

with open("ligand_network.graphml", mode='w') as f:
    f.write(ligand_network.to_graphml())

# 2.3 Create the relationship of the mapping 

output_dir = Path("2d_transformations")
output_dir.mkdir(exist_ok=True)

print("Saving figures transformations ...")

network_edges = {
    f"{edge.componentA.name} ➔ {edge.componentB.name}": edge 
    for edge in ligand_network.edges
}

for edge in ligand_network.edges:
    # Names of (SVG o PNG)
    filename = output_dir / f"mapping_{edge.componentA.name}_to_{edge.componentB.name}.svg"

    # Save Figure
    edge.draw_to_file(filename)

print(f"\n All transformations images saved in {output_dir}")


# In[104]:


# 2.4 Figure with scores of the network
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
def export_network_svg(net, filename="ligand_network_scores.svg"):
    G = nx.Graph()
    for m in net.edges:
        s = m.annotations.get("score")
        G.add_edge(
            m.componentA.name,
            m.componentB.name,
            score=f"{float(s):.2f}" if s else "",
        )

    # Layout 
    pos = nx.kamada_kawai_layout(G)
    nodes = list(pos.keys())
    for _ in range(50):
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                p1, p2 = np.array(pos[nodes[i]]), np.array(pos[nodes[j]])
                d = np.linalg.norm(p1 - p2)
                if 0 < d < 0.22:
                    shift = ((0.22 - d) / 2) * (p1 - p2) / d
                    pos[nodes[i]] += shift
                    pos[nodes[j]] -= shift
    fig, ax = plt.subplots(figsize=(12, 8))

    nx.draw_networkx_edges(
        G, pos, edge_color="#555555", width=1.3, alpha=0.9, ax=ax
    )
    nx.draw_networkx_labels(
        G,
        pos,
        font_size=10,
        font_weight="bold",
        font_color="#1F4E79",
        bbox=dict(boxstyle="round,pad=0.3", fc="w", ec="none"),
        ax=ax,
    )
    nx.draw_networkx_edge_labels(
        G,
        pos,
        edge_labels=nx.get_edge_attributes(G, "score"),
        font_size=9,
        rotate=False,
        bbox=dict(boxstyle="round,pad=0.2", fc="w", ec="#D0D0D0"),
        ax=ax,
    )

    plt.title(
        "Ligand Network with Scores", fontsize=14, pad=15, fontweight="bold"
    )
    plt.axis("off")
    plt.tight_layout()


    plt.savefig(filename, bbox_inches="tight")
    plt.close()

export_network_svg(ligand_network)


# In[105]:


# 3.1 Define the components of the system

solvent = openfe.SolventComponent() #Default 0.15 M
protein = openfe.ProteinComponent.from_pdb_file(pdb_prot)

mapping = next(iter(ligand_network.edges))

# The molecules of the transformation
systemA = openfe.ChemicalSystem({
    'ligand': mapping.componentA,
    'solvent': solvent,
    'protein': protein
})
systemB = openfe.ChemicalSystem({
    'ligand': mapping.componentB,
    'solvent': solvent,
    'protein': protein    
})


# In[ ]:


# 4.1 Settings of the protocol 
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
settings = RelativeHybridTopologyProtocol.default_settings() # This loads the default settings for a protocol

# 4.2 Change settings (If you are not aware of the things you are changing let the default settings )
from openff.units import unit

# 4.3 Special characteristics of solvent and complex systems 
#------------------------------------------------------------
#               ADVANCED SETTINGS (BASE)
#------------------------------------------------------------
##Temperature
settings.thermo_settings.temperature = 310.0 * unit.kelvin #Temperature
##Forcefields (Theres a bunch more of FF, but amber14 is the one tested in the paper)
#settings.forcefield_settings.forcefields = [   
#    'amber/ff14SB.xml', # Options amber/ff14SB.xml
#    'amber/tip3p_standard.xml', # Options amber19/tip3p_standard.xml,  amber19/tip4pew.xml
#    'amber/tip3p_HFE_multivalent.xml', # Remove comment if using metals with tip3p
#    'amber/phosaa10.xml'  # For PTMs 
#]
##Water model
#settings.solvation_settings.solvent_model = 'tip3p' #opc model is not included
##Box_shape
#settings.solvation_settings.box_shape = 'dodecahedron' # Options cube, dodecahedron
##FF small molecule
settings.forcefield_settings.small_molecule_forcefield = 'openff-2.2.1'
##Equilibration
#settings.simulation_settings.equilibration_length = 1.0 * unit.nanosecond    # EQUILIBRATION 
##Simulations 
#settings.simulation_settings.production_length = 5.0 * unit.nanosecond     #PRODUCTION
##Lambda windows
#settings.simulation_settings.n_replicas = 11      # 
#settings.lambda_settings.lambda_windows = 11      # replicas and lambda must be the same
##Replicas
#settings.protocol_repeats = 3
#settings.simulation_settings.sampler_method = "repex" # Sampler method, options are: repex (H-REX), sams (SAMS), independent (each window a simulations)
#------------------------------------------------------------
#       COPY OF BASE SETTINGS AND SPECIAL PARAMETERS
#------------------------------------------------------------
#Solvent
solvent_settings = settings.model_copy(deep=True)
#Complex
complex_settings = settings.model_copy(deep=True)
complex_settings.solvation_settings.solvent_padding = 1 * unit.nanometer

#------------------------------------------------------------
#                   SYSTEMS SETTINGS
#------------------------------------------------------------
## Solvent
solvent_protocol = RelativeHybridTopologyProtocol(solvent_settings)
## Complex
complex_protocol = RelativeHybridTopologyProtocol(complex_settings)



# In[147]:


##See an special setting
#complex_settings
#solvent_settings
complex_settings.simulation_settings.sampler_method


# In[124]:


# 4.4 Save the protocols if you want to make sure hows going to run
import yaml
from pathlib import Path

def settings_to_dict(settings_obj):
    """Used Pydantic to convert"""
    if hasattr(settings_obj, "model_dump"):
        return settings_obj.model_dump()
    return settings_obj.dict()

def export_yaml(settings_obj, filename):
    filepath = Path("protocol_configs") / filename
    filepath.parent.mkdir(exist_ok=True)

    # Get dictionary
    data = settings_to_dict(settings_obj)

    with open(filepath, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    print(f"Protocols saved in: {filepath}")

# Export
export_yaml(complex_settings, "complex_protocol.yaml")
export_yaml(solvent_settings, "solvent_protocol.yaml")


# In[125]:


# 4.4 Definition of all transformations
transformations = []
for mapping in ligand_network.edges:
    for leg in ['solvent', 'complex']:
        # use the solvent and protein created above
        sysA_dict = {'ligand': mapping.componentA,
                     'solvent': solvent}
        sysB_dict = {'ligand': mapping.componentB,
                     'solvent': solvent}

        if leg == 'complex':
            # If this is a complex transformation we use the complex protocol
            # and add in the protein to the chemical states
            protocol = complex_protocol
            sysA_dict['protein'] = protein
            sysB_dict['protein'] = protein
        else:
            # If this is a solvent transformation we just use the solvent protocol
            protocol = solvent_protocol

        # we don't have to name objects, but it can make things (like filenames) more convenient
        sysA = openfe.ChemicalSystem(sysA_dict, name=f"{mapping.componentA.name}_{leg}")
        sysB = openfe.ChemicalSystem(sysB_dict, name=f"{mapping.componentB.name}_{leg}")

        prefix = "rbfe_"  # prefix is only to exactly reproduce CLI

        transformation = openfe.Transformation(
            stateA=sysA,
            stateB=sysB,
            mapping=mapping,
            protocol=protocol,  # use protocol created above
            name=f"{prefix}{sysA.name}_{sysB.name}"
        )
        transformations.append(transformation)

network = openfe.AlchemicalNetwork(transformations)

import pathlib
# First we create the directory
transformation_dir = pathlib.Path("transformations")
transformation_dir.mkdir(exist_ok=True)

# then we write out each transformation
for transformation in network.edges:
    transformation.to_json(transformation_dir / f"{transformation.name}.json")


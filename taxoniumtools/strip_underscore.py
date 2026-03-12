import dendropy
import os

input_file = '/home/vlad/code/taxonium/taxoniumtools/src/taxoniumtools/testtree/HA_tree_rerooted.nwk'
output_file = '/home/vlad/code/taxonium/taxoniumtools/src/taxoniumtools/testtree/HA_tree_stripped.nwk'

# Load the tree. preserve_underscores=True is often useful in Newick
tree = dendropy.Tree.get(path=input_file,
                         schema='newick',
                         preserve_underscores=True)

# Iterate over all taxa and rename them
for taxon in tree.taxon_namespace:
    if taxon.label and '_' in taxon.label:
        parts = taxon.label.rsplit('_', 1)
        taxon.label = parts[0]
        print(f"Renamed: {taxon.label}_{parts[1]} -> {taxon.label}")

# Also check node labels for internal nodes if they are named
for node in tree.nodes():
    if node.label and '_' in node.label:
        parts = node.label.rsplit('_', 1)
        node.label = parts[0]
        print(f"Renamed node: {node.label}_{parts[1]} -> {node.label}")

# Save the modified tree
tree.write(path=output_file, schema='newick')
print(f"Saved to {output_file}")

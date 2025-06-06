import ROOT
import argparse

input_file = "/Users/yfeng/Desktop/TTU/CaloX/Data/run316_250517140056.root"
output_file = "/Users/yfeng/Desktop/TTU/CaloX/Data/run316_250517140056_filtered.root"
# Open the original ROOT file
input_file = ROOT.TFile(input_file, "READ")

# Get the tree from the input file
input_tree = input_file.Get("EventTree")

# Create a new ROOT file to store the filtered tree
output_file = ROOT.TFile(output_file, "RECREATE")

# Create a new tree in the output file with the same structure as the input tree
output_tree = input_tree.CloneTree(0)

# Loop over the entries in the input tree and copy the desired ones
for entry in input_tree:
    if entry.FERS_Board5_energyHG[0] > 2000:
        output_tree.Fill()

# Write the new tree to the output file
output_tree.Write()

# Close both files
input_file.Close()
output_file.Close()

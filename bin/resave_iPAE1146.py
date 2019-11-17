import cobra

# Load the PAO1 model
pa = cobra.io.read_sbml_model('../data/previous_reconstructions/iPAE1146.xml')

# save the PAO1 model, which should no longer raise thousands of SBML errors
cobra.io.write_sbml_model(pa,'../data/previous_reconstructions/iPAE1146_resaved.xml')
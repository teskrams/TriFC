from TriFC.IfcFile import IfcFile

ifcfile = IfcFile("/Users/theodor/kode/master/TriFC/triFC/example_ifc/academic_test_2.ifc")

# #### Setting the parameters if you are not running in academic mode ########################
# ############################################################################################
# ifcfile.set_ifc_file_read_config(
#     pset_quality = "01 Merknader",                      # Property set for material quality
#     quality_attribute = "K01 Materialkvalitet",         # Atteribute for material quality
#     pset_prod_method = "01 Merknader",                  # Property set for production method
#     prod_method_attribute = "K11 Produksjonsmetode",    # Atteribute for production method
#     degree = 60                                         # Maximum slope Formwork
#     )
# ##############################################################################################

ifcfile.run_full_calculation(
    name="grubb.csv", 
    result_filepath="/Users/theodor/kode/master/TriFC/triFC/example_results",  
    ifc_name="academic_test_2_TriFC.ifc", 
    ifc_path="/Users/theodor/kode/master/TriFC/triFC/example_ifc", 
    printing=2, 
    academic=True
    )

#### Example of how to use the visualization tool on the calulated elements ################
############################################################################################
# from TriFC.IfcPlotter import IfcPlotter

# ifc_plotter = IfcPlotter()

# ifc_plotter.interactive_plot(
#     elements=[ifcfile.elements[1]],                     # List of IfcElements
#     collisions=[ifcfile.elements[1].formwork],          # Bool lists ref tri in elemnts
#     #check=1,                                            # Coloring tri at id
#     padding=0.5,                                        # Padding around element
#     random_color=False                                  # Random color triangles
#     )
##############################################################################################
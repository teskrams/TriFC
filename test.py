from Ifc import IfcFile

ifcfile = IfcFile("/Users/theodor/kode/master/TriFC/triFC/example_ifc/academic_test_2.ifc")

ifcfile.run_full_calculation(name="academic_test_2.csv", 
                             result_filepath="/Users/theodor/kode/master/TriFC/triFC/example_results",  
                             ifc_name="academic_test_2_runned_through_model.ifc", 
                             ifc_path="/Users/theodor/kode/master/TriFC/triFC/example_ifc", 
                             printing=2, 
                             academic=True
                             )
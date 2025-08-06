import pihmmi
from config import (
    prot_dir,
    df, df_class_1, df_class_2, df_class_3, df_class_4,df_class_2v2,
    dataset123, dataset1, dataset2, dataset3, dataset4,
    gen_name_phac,file_name_phac,search_file,
    uniport_url_phac4,orthodb_url_phac4,uniport_url_file_name,orthodb_url_file_name
)



def main(dataset=False,num=1):
    if dataset:
        pihmmi.create_dataset(output_dir=prot_dir,
                              dataset_name=dataset123,
                              gen_name_list=gen_name_phac,
                              file_name_list=file_name_phac,
                              uniport_url=uniport_url_phac4,
                              orthodb_url=orthodb_url_phac4,
                              additional_uniport_file=uniport_url_file_name,
                              additional_orthodb_file=orthodb_url_file_name)

    df_classes = [df_class_1, df_class_2, df_class_3, df_class_4]

    # PhaC I
    if num==1:
        results = pihmmi.run_pihmmi_pipeline(
            df_class=df_classes[0],
            number_pipeline=1,
            email="247034@vutbr.cz",
            model_name=search_file,
            iterations=3,
            max_seq=10,
            multi_alignment=True,
            silhouette_analysis=False,
            data_set=dataset4,
            output_dir=prot_dir
        )

    # PhaC II
    if num==2:
        results = pihmmi.run_pihmmi_pipeline(
            df_class=df_class_2,
            number_pipeline=2,
            email="247034@vutbr.cz",
            model_name=search_file,
            iterations=3,
            max_seq=10,
            multi_alignment=True,
            cluster_tree=True,
            silhouette_analysis=False,
            data_set=dataset4,
            output_dir=prot_dir
        )
    # PhaC II v2
    if num==3:
        results = pihmmi.run_pihmmi_pipeline(
            df_class=df_class_2v2,
            number_pipeline=2,
            email="247034@vutbr.cz",
            model_name=search_file,
            iterations=3,
            max_seq=10,
            multi_alignment=True,
            cluster_tree=True,
            silhouette_analysis=False,
            data_set=dataset4,
            output_dir=prot_dir
        )
    # PhaC III
    if num==4:
        results = pihmmi.run_pihmmi_pipeline(
            df_class=df_classes[2],
            number_pipeline=3,
            email="247034@vutbr.cz",
            model_name=search_file,
            iterations=5,
            max_seq=7,
            multi_alignment=False,
            cluster_tree=False,
            silhouette_analysis=True,
            data_set=dataset4,
            output_dir=prot_dir
        )

    # PhaC IV
    if num == 5:
        results = pihmmi.run_pihmmi_pipeline(
            df_class=df_class_4,
            number_pipeline=4,
            email="247034@vutbr.cz",
            model_name=search_file,
            iterations=5,
            max_seq=7,
            multi_alignment=False,
            cluster_tree=True,
            silhouette_analysis=True,
            data_set=dataset4,
            output_dir=prot_dir
        )

    return results

if __name__ == "__main__":
    # #PhaC I
    # main(dataset=True,num=1)
    # # PhaC II
    # main(dataset=False, num=2)
    #PhaC II v2
    main(dataset=False,num=3)
    # #PhaC III
    # main(dataset=False,num=4)
    # #PhaC IV
    # main(dataset=False,num=5)





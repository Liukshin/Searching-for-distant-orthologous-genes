import pihmmi
from config import (
    prot_dir,
    df, df_class_1, df_class_2, df_class_3, df_class_4,df_class_2v2,
    dataset123, dataset1, dataset2, dataset3, dataset4,
    gen_name_phac,file_name_phac,search_file,
    uniport_url_phac4,orthodb_url_phac4,uniport_url_file_name,orthodb_url_file_name
)




pihmmi.create_dataset(output_dir=prot_dir,
                      dataset_name=dataset123,
                      gen_name_list=gen_name_phac,
                      file_name_list=file_name_phac,
                      uniport_url=uniport_url_phac4,
                      orthodb_url=orthodb_url_phac4,
                      additional_uniport_file=uniport_url_file_name,
                      additional_orthodb_file=orthodb_url_file_name)


df_classes = [df_class_1, df_class_2, df_class_3, df_class_4]
#PhaC I
# results = pihmmi.run_pihmmi_pipeline(
#     df_class=df_classes[0],
#     number_pipeline=1,
#     email="247034@vutbr.cz",
#     iterations=3,
#     max_seq=10,
#     multi_alignment=True,
#     silhouette_analysis = True,
#     data_set=dataset4,
#     output_dir=prot_dir
# )

#PhaC III

# results = pihmmi.run_pihmmi_pipeline(
#     df_class=df_classes[2],
#     number_pipeline=3,
#     email="247034@vutbr.cz",
#     iterations=5,
#     max_seq=7,
#     multi_alignment=False,
#     cluster_tree= False,
#     silhouette_analysis = True,
#     data_set=dataset4,
#     output_dir=prot_dir
# )

#phac
results = pihmmi.run_pihmmi_pipeline(
    df_class=df_class_4,
    number_pipeline=4,
    email="247034@vutbr.cz",
    model_name=search_file,
    iterations=5,
    max_seq=7,
    multi_alignment=False,
    cluster_tree= True,
    silhouette_analysis = True,
    data_set=dataset4,
    output_dir=prot_dir
)


print("Pipeline Results:")
print(f"- Found sequences: {len(results.get('found_sequences', []))}")
print(f"- Clusters: {len(results.get('clusters', {}))}")


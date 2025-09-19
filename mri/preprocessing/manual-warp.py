import os
import aux_warp
import subprocess

def apply_transform_segmentation(input_img, output_img, ref_img, transformation):
    # input_img: segmentation file from atlas
    # output_img: new segmentation specified with animal_id
    # ref_img: unregistered, resampled TurboRARE img of animal_id
    # transformation: transformation matrix output of ANTs registration

    s1 = ['antsApplyTransforms', '--default-value', '0', '--float', '0', '--input']
    s2 = ['--input-image-type', '0', '--interpolation', 'BSpline[5]', '--output']
    s3 = ['--reference-image']
    s4 = ['--transform']
    s5 = ["--interpolation", "GenericLabel"]

    cmnd = s1 + [input_img] + s2 + [output_img] + s3 + [ref_img] + s4 + ["[" + transformation + ", 1]"] + s5
    result = subprocess.run(cmnd, capture_output=True, text=True)
    print(result.stdout)



root="/Users/eminhanozil/Dropbox (Yanik Lab)/Localization Manuscript 2024/RAT DATA"

dwi_path = os.path.join(root, "WHS_SD_rat_atlas_v4_pack","WHS_SD_rat_DWI_v1.01.nii.gz")
mask_path = os.path.join(root, "WHS_SD_rat_atlas_v4_pack","WHS_SD_v2_brainmask_bin.nii.gz")
t2s_path = os.path.join(root, "WHS_SD_rat_atlas_v4_pack","WHS_SD_rat_T2star_v1.01.nii.gz")
atlas_path = os.path.join(root, "WHS_SD_rat_atlas_v4_pack","WHS_SD_rat_atlas_v4.nii.gz")

animal_id="rEO_10"
session="postsurgery20241024"
bids_work="_ind_type_40"

filename_segmentation = "relaxation_heatmap_resampled.nii.gz"

input_segmentation = os.path.join(root, animal_id, "mri", "ses-"+session, "anat", filename_segmentation)
transformation = os.path.join(root, animal_id, "mri", "ses-"+session, bids_work, "s_register", "output_Composite.h5")

output_filename = "heatmap_registered.nii.gz"
output_segmentation = os.path.join(root, animal_id, "mri", "ses-"+session, "anat", output_filename)

result = apply_transform_segmentation(input_segmentation, output_segmentation, t2s_path, transformation)





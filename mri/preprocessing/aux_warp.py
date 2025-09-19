import subprocess

def biascorrect(input_img, output_img):
    s1 = ['N4BiasFieldCorrection', '--bspline-fitting', '[10,4]', '-d', '3', '--input-image']
    s2 = ['--convergence', '[150x100x50x30,1e-16]', '--output']
    s3 = ['--shrink-factor', '2']

    cmnd = s1 + [input_img] + s2 + [output_img] + s3
    result = subprocess.run(cmnd, capture_output=True, text=True)
    print(result.stdout)


def apply_transform(input_img, output_img, ref_img, transformation):
    # input_img: animal high_res or axial postsurgery img, pre-registered manually to presurgery space
    # output_img: registered high_res or axial postsurgery img to atlas space
    # ref_img: atlas

    s1 = ['antsApplyTransforms', '--default-value', '0', '--float', '0', '--input']
    s2 = ['--input-image-type', '0', '--interpolation', 'BSpline[5]', '--output']
    s3 = ['--reference-image']
    s4 = ['--transform']

    cmnd = s1 + [input_img] + s2 + [output_img] + s3 + [ref_img] + s4 + [transformation]
    result = subprocess.run(cmnd, capture_output=True, text=True)
    print(result.stdout)


def apply_transform_atlas(input_img, output_img, ref_img, transformation):
    # input_img: atlas
    # output_img: new segmentation specified with animal_id
    # ref_img: unregistered, resampled TurboRARE img of animal_id

    s1 = ['antsApplyTransforms', '--default-value', '0', '--float', '0', '--input']
    s2 = ['--input-image-type', '0', '--interpolation', 'BSpline[5]', '--output']
    s3 = ['--reference-image']
    s4 = ['--transform']
    #     s5=["--interpolation", "GenericLabel"]

    cmnd = s1 + [input_img] + s2 + [output_img] + s3 + [ref_img] + s4 + ["[" + transformation + ", 1]"]
    result = subprocess.run(cmnd, capture_output=True, text=True)
    print(cmnd)
    print(result.stderr)
    print(result.stdout)
    return result


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

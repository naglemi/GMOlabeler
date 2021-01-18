from PIL import Image
import glob
import os
import sys

def get_concat_h(im1, im2): #https://note.nkmk.me/en/python-pillow-concat-images/
    dst = Image.new('RGB', (im1.width + im2.width, im1.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (im1.width, 0))
    return dst
    
def main(wd, alpha, mode, rotate, angle, chroma_standard_prefix = 'hroma', verbose = True):
    os.chdir(wd)
    print(sys.argv)
    #directory_contents = os.listdir('.')
    segmentation_output_paths = glob.glob('*segment_uncropped.png')
    rgb_paths = glob.glob('*rgb.jpg')
    rgb_paths = [i for i in rgb_paths if not (str(chroma_standard_prefix) in i)]
    rgb_paths.sort()
    segmentation_output_paths.sort()
    if len(rgb_paths) != len(segmentation_output_paths):
        print('Error: There is not a segmentation output for every rgb image.')
        print(len(rgb_paths))
        print(len(segmentation_output_paths))
        sys.exit(1)
    for i in range(1, len(rgb_paths)):
        if verbose == True:
            print('Building composite for image ' + rgb_paths[i])
        segment_img = Image.open(segmentation_output_paths[i]).convert("RGBA").rotate(180)
        rgb_img = Image.open(rgb_paths[i]).convert("RGBA").rotate(180)
    
        if rotate == True:
            rgb_img = rgb_img.rotate(angle)
            segment_img = segment_img.rotate(angle)
        
        image_blended = Image.blend(segment_img, rgb_img, alpha=float(alpha))

        if (mode == 'composite' or mode == 'both'):
            composite = get_concat_h(rgb_img, segment_img)
            composite = get_concat_h(composite, image_blended)
            output_filename = rgb_paths[i].replace('_rgb.jpg', '_composite.png')
            composite.save(output_filename)
            #composite.show()
        if (mode == 'blend' or mode == 'both'):
            output_filename = rgb_paths[i].replace('_rgb.jpg', '_blend.png')
            image_blended.save(output_filename)
            ##composite.show()
        
if __name__== "__main__":
    main(wd = sys.argv[1],
    alpha = sys.argv[2], 
    mode = sys.argv[3],
    rotate = sys.argv[4],
    angle = sys.argv[5])

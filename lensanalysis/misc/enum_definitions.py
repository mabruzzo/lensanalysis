from enum import Enum

# in python 3 it would be better to use the Flag subclass of enum
class DescriptorEnum(Enum):
    tomo = 1
    smoothed = 2
    noisy = 3

# I am not defining analysis map options as a tuple because this list will
# in principle be added to 
analysis_object_names = ["conv_map",
                         "shear_map",
                         "shear_cat",
                         "peak_loc",
                         "peak_counts"]

# gives the valid descriptors for the different analysis objects
analysis_object_descriptors = {"conv_map" : (DescriptorEnum.tomo,
                                             DescriptorEnum.smoothed,
                                             DescriptorEnum.noisy),
                               "shear_map" : (DescriptorEnum.tomo,
                                              DescriptorEnum.noisy),
                               "shear_cat" : (),
                               "peak_loc" : (DescriptorEnum.tomo,),
                               "peak_counts" : (DescriptorEnum.tomo,)}

descriptor_mapping = {DescriptorEnum.tomo : ["tomo","tomographic"],
                      DescriptorEnum.smoothed : ["smooth","smoothed"],
                      DescriptorEnum.noisy : ["noisy"]}

from itertools import chain, combinations, izip_longest

from enum import Enum



# in python 3 it would be better to use the Flag subclass of enum
class DescriptorEnum(Enum):
    tomo = 1
    smoothed = 2
    noisy = 3

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

def all_combinations(omit_analysis_objects = [], omit_descriptors = []):
    out = ()

    for key,value in analysis_object_descriptors.iteritems():
        if key not in omit_analysis_objects:
            cur = (((),key),)
            descriptors = [descr for descr in value if descr not in
                           omit_descriptors]
            #for descr in value:
            #    if descr in omit_descriptors:
            #        continue
            #    else:
            #        descriptors.append(descr)
            for i in range(1,len(descriptors)+1):
                temp = combinations(descriptors,i)
                temp = izip_longest(temp, (key,), fillvalue = key)
                cur = chain(cur,temp)
        out = chain(out,cur)
    return out

if __name__ == '__main__':
    print list(all_combinations(["shear_cat"]))


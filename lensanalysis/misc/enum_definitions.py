from itertools import chain, combinations, izip_longest

from aenum import Flag

class Descriptor(Flag):
    none = 0
    tomo = 1
    smoothed = 2
    noisy = 4
    tomo_smoothed = tomo | smoothed
    tomo_noisy = tomo | noisy
    smoothed_noisy = smoothed | noisy
    tomo_smoothed_noisy = tomo | smoothed_noisy

# below are the descriptors of the 3 fundamental flags
single_descriptor_flags = [("tomo", Descriptor.tomo),
                            ("smoothed", Descriptor.smoothed),
                            ("noisy",Descriptor.noisy)]

# gives the valid descriptors for the different analysis objects
analysis_object_descriptors = {"conv_map" : (Descriptor.tomo,
                                             Descriptor.smoothed,
                                             Descriptor.noisy),
                               "shear_map" : (Descriptor.tomo,
                                              Descriptor.noisy),
                               "shear_cat" : (),
                               "peak_loc" : (Descriptor.tomo,),
                               "peak_counts" : (Descriptor.tomo,)}

descriptor_mapping = {Descriptor.tomo : ["tomo","tomographic"],
                      Descriptor.smoothed : ["smooth","smoothed"],
                      Descriptor.noisy : ["noisy"],
                      Descriptor.none : []}

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

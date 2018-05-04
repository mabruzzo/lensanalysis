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
                                             Descriptor.noisy,
                                             Descriptor.none),
                               "shear_map" : (Descriptor.tomo,
                                              Descriptor.noisy,
                                              Descriptor.none),
                               "shear_cat" : (Descriptor.none),
                               "peak_loc" : (Descriptor.tomo,
                                             Descriptor.none),
                               "peak_counts" : (Descriptor.tomo,
                                                Descriptor.none),
                               "power_spectrum" : (Descriptor.tomo)}

descriptor_mapping = {Descriptor.tomo : ["tomo","tomographic"],
                      Descriptor.smoothed : ["smooth","smoothed"],
                      Descriptor.noisy : ["noisy"],
                      Descriptor.none : []}


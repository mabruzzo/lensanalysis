from procedure import AtomicProcedureStep, ConversionProcedureStep


class LoadCollection(ConversionProcedureStep):
    """
    Loads the collection entries.

    In effect, we are converting an integer representing the realization to an 
    instance of the collection.
    """

    def __init__(self,collection_storage):
        self.collection_storage

class SaveCollectionEntries(AtomicProcedureStep):
    """
    Configurable step that saves the entries of a collections.
    """

    def __init__(self,saver,fname_formatter,attribute_name):
        self._saver = saver
        self._fname_formatter = fname_formatter
        self._attribute_name = attribute_name

    def apply(self,data_object, packet):
        """
        Probably will need to modify once we implement some of the actual 
        objects.
        """
        data_id = packet.data_id

        for entry,fname in self._fname_formatter.format_pairs(data_object,
                                                              data_id):
            saver.save(entry,fname)

from procedure import AtomicProcedureStep, ConversionProcedureStep
from ..misc.serialization import CollectionSaver, CollectionLoader
from ..misc.log import logprocedure

class LoadCollection(ConversionProcedureStep):
    """
    Loads the collection entries.

    In effect, we are converting an integer representing the realization to an 
    instance of the collection.

    Parameters
    ----------
    collection_loader : CollectionLoader
        An instance of a (virtual) subclass of abstract base class 
        CollectionLoader.
    """

    def __init__(self,collection_loader):
        if isinstance(collection_loader, CollectionLoader):
            self.collection_loader = collection_loader
        else:
            raise TypeError("collection_storage must be an instance of a "
                            "(virtual) subclass of abc, CollectionStorage")

    def conversion_operation(self,data_object,packet):
        logprocedure.debug("Loading realization {:d}".format(packet.data_id))
        return self.collection_loader.load(packet.data_id)

class SaveCollectionEntries(AtomicProcedureStep):
    """
    Configurable step that saves the entries of a collections.

    Parameters
    ----------
    collection_saver : CollectionSaver
        An instance of a (virtual) subclass of abstract base class 
        CollectionSaver.
    """

    def __init__(self,collection_saver):
        if isinstance(collection_saver, CollectionSaver):
            self.collection_saver = collection_saver
        else:
            raise TypeError("collection_storage must be an instance of a "
                            "(virtual) subclass of abc, CollectionStorage")

    def apply(self,data_object, packet):
        """
        Probably will need to modify once we implement some of the actual 
        objects.
        """
        logprocedure.debug("Saving realization {:d}".format(packet.data_id))
        self.collection_saver.save(data_object,packet.data_id)

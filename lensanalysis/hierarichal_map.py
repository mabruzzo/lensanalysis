from collections import Sequence, MutableMapping
import copy
import itertools
    

class HierarichalMapping(MutableMapping):
    """
    Tracks hierarichal data.

    Because the map is mutable, we do not treat slices as views
    """

    def __init__(self,fields):
        if isinstance(fields,basestring):
            fields = [fields]
        elif not isinstance(fields,Sequence):
            raise ValueError("fields must be a str or a sequence of str")
        self._fields = fields
        self._mapping = {}

    @property
    def fields(self):
        return self._fields

    def keys(self):
        return sorted(self._mapping.keys())

    def as_dict(self):
        # I believe we actually want a shallow copy
        return copy.copy(self._mapping)
    
    def __getitem__(self,key):
        """
        The entries are taken in order of the 
        """
        if not isinstance(key,tuple):
            key = (key,)

        for elem in key:
            if any([type(elem)==type(slice) for elem in key]):
                raise NotImplementedError("Cannot yet handle slicing")
            
        if len(key) == len(self._fields):
            return self._mapping[key]
        elif len(key) < len(self._fields):
            temp = list(key)
            for elem in range(len(self._fields)-len(key)):
                temp.append(slice)
            key = tuple(temp)
        else:
            raise ValueError("length of key must be <= number of fields")

        remaining_entries = self.keys()

        slice_loc = [False for elem in key]
        for i,elem in enumerate(key):
            if type(elem)==type(slice):
                slice_loc[i] = True
            else:
                temp = []
                for x in remaining_entries:
                    print elem
                    if x[i] == elem:
                        temp.append(x)
                remaining_entries = temp


        new_fields = [field for field,keep in zip(self.fields,
                                                  slice_loc) if keep]

        new_mapping = HierarichalMapping(new_fields)
        for elem in remaining_entries:
            new_key = tuple([val for val,keep in zip(elem,
                                                     slice_loc) if keep])
            new_mapping[new_key] = self[elem]

        return new_mapping

    def __setitem__(self, key, item):
        if not isinstance(key,tuple):
            key = (key,)
        elif len(key)==len(self._fields):
            self._mapping[key]=item
        else:
            raise ValueError("length of key must be equal to number of fields")

    def __delitem__(self,key,item):
        if not isinstance(key,tuple):
            key = (key,)
        elif len(key) == len(self._fields):
            del self._mapping[key]
        else:
            raise ValueError("length of key must be equal to number of fields")

    def __iter__(self):
        for elem in self.keys():
            yield elem

    def __len__(self):
        return len(self._mapping.keys())


if __name__ == '__main__':

    t = {"f" : 5,
         2 : 4}
    
    
    t = HierarichalMapping(["field0","field1","field2"])
    t["a",3,4] =5
    t["a",7,9] = 10
    b = t["a"]
    print b.keys()
    print b.fields
    for elem in b:
        print elem,b[elem]

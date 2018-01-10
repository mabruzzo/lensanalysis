from abc import ABCMeta
from collections import Iterable

import numpy as np

class DataPacket(object):
    def __init__(self,data_id,features = []):
        self.data_id = data_id
        self.features = np.array(features)
        

class ComponentProceudureStep(object):
    __metaclass__ = ABCMeta
    """
    Abstract base class to represent a procedural step

    Procedures will be implemented using the Composite design pattern
    """

    @abstractmethod
    def apply(self, data_object, packet):
        """
        Apply the procedural step to the packet.
        """
        pass

    @classmethod
    def __subclasshook__(cls,C):
        if cls is ComponentProceudureStep:
            if any("apply"  in B.__dict__ for B in C.__mro__):
                return True
        return NotImplemented


class CompositeProcedureStep(ComponentProcedureStep):
    """
    Procedural Step constructed from other Procedure Steps
    """

    def __init__(self):
        self._child_procedure_steps = []

    def add(self, procedure_step):
        if not isinstance(procedural_step, ComponentProceudureStep):
            raise ValueError("Only accepts obects with apply method")
        self._child_procedure_steps.append(procedure_step)

    def remove(self, procedure_step):
        self._child_procedure_steps.remove(procedure_step)

    def apply(self, data_object, packet):
        for procedure_step in self._child_procedure_steps:
            proedure_step.apply(data_object, packet)

class Procedure(CompositeProcedureStep):
    __metaclass__ = ABCMeta
    """
    Base class to represent the procedure

    It acts exactly like an instance of CompositeProcedureStep plus it can
    be used to apply a procedure to an iterable.

    Notes:
    Could subclass to provide a CheckpointedProcedure, that allows you to skip 
    over previously checkpointed steps/
    """

    def apply_to_iterable(self,iterable):
        if not isinstance(iterable, Iterable):
            raise ValueError("iterable must be an instance of a class with the "
                             "methods of the Iterable Abstract Base Class")
        for data_object, packet in iterable:
            self.apply(data_object, packet)

class SimpleProcedure(Procedure):
    """
    Simplest concrete implementation of Procedure
    """

class ProcedureStepDecorator(ComponentProcedureStep):
    __metaclass__ = ABCMeta

    def get_wrapped_step(self):
        pass

    def set_wrapped_step(self,value):
        pass

    wrapped_step = abstractproperty(get_wrapped_step,set_wrapped_step)

class IntermediateProcedureStep(ProcedureStepDecorator):
    __metaclass__ = ABCMeta
    """
    A procedure step where the output can be fed into the next step
    (Example: smoothing a convergence Map OR loading an object
    This is very similar to AbstractConversionProcedureStep - except the 
    abstract type of the data object is not changing. If I choose to keep both 
    classes, then AbstractConversionProcedureStep will be the subclass of this 
    class.

    This would wrap around another procedure step.
    """


    @abstractmethod
    def intermediate_operation(self,data_object,packet):
        pass

    def apply(self, data_object,packet):
        result = self.intermediate_operation(data_object,packet)
        self._apply_to_results(result,packet)
        if wrapped_step is not None:
            self.wrapped_step(result,packet)


"""
The building of a procedure is going to ordinarily revolve around Conversion 
Procedure Steps.
"""

class ConversionProcedureStep(ProcedureStepDecorator):
    __metaclass__ = ABCMeta

    """
    Handles operations that convert data objects between different abstract 
    types.

    The conversion operation has been factored out to the abstract 
    conversion_operation and the application of further steps has been 
    factored out to the abstract wrapped_step property
    """

    @abstractmethod
    def conversion_operation(self,data_object,packet):
        pass

    def apply(self, data_object,packet):
        result = self.conversion_operation(data_object,packet)
        self._apply_to_results(result,packet)
        if wrapped_step is not None:
            self.wrapped_step(result,packet)

    @classmethod
    def __subclasshook__(cls,C):
        raise NotImplementedError

    
class ConversionProcedureStepGrouping(object):
    """
    Grouping of steps involving the result of the conversion of data objects 
    between different abstract types
    
    Probably not the best way to go about implementing this idea, but the idea 
    is that after converting between abstract types there are 2 classes of 
    operations applied to the result:
    1. processing of the product - saving, computing features
    2. conversion to other abstract data types

    The idea would be for this object to wrap around an instance of an 
    implementation of a ConversionProcedureStep object, and manage the 
    procedure step wrapped by that instance.
    Specifically, that instance would wrap a CompositeProcedureStep object and 
    an instance of this class would manage the entries in that object
    """

        
class AtomicProcedureStep(ComponentProcedureStep):
    """
    Base class for a procedure step that is not made up of other steps

    This is the Leaf in the composite design
    """

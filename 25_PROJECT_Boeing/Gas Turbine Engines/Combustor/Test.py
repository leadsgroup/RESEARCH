import matplotlib.pyplot as plt
import yaml

   
class Solution(Transport, Kinetics, ThermoPhase):
    """
    A class for chemically-reacting solutions. Instances can be created to
    represent any type of solution -- a mixture of gases, a liquid solution, or
    a solid solution, for example.

    Class `Solution` derives from classes `ThermoPhase`, `Kinetics`, and
    `Transport`.  It defines no methods of its own, and is provided so that a
    single object can be used to compute thermodynamic, kinetic, and transport
    properties of a solution.

    To skip initialization of the Transport object, pass the keyword argument
    ``transport_model=None`` to the `Solution` constructor.

    The most common way to instantiate `Solution` objects is by using a phase
    definition, species and reactions defined in an input file::

        gas = ct.Solution('gri30.yaml')

    If an input file defines multiple phases, the corresponding key in the
    ``phases`` map can be used to specify the desired phase via the ``name`` keyword
    argument of the constructor::

        gas = ct.Solution('diamond.yaml', name='gas')
        diamond = ct.Solution('diamond.yaml', name='diamond')

    The name of the `Solution` object defaults to the *phase* identifier
    specified in the input file. Upon initialization of a `Solution` object,
    a custom name can assigned via::

        gas.name = 'my_custom_name'

    `Solution` objects can also be constructed using `Species` and `Reaction`
    objects which can themselves either be imported from input files or defined
    directly in Python::

        spec = ct.Species.list_from_file("gri30.yaml")
        spec_gas = ct.Solution(thermo='ideal-gas', species=spec)
        rxns = ct.Reaction.list_from_file("gri30.yaml", spec_gas)
        gas = ct.Solution(thermo='ideal-tas', kinetics='gas',
                          species=spec, reactions=rxns, name='my_custom_name')

    where the ``thermo`` and ``kinetics`` keyword arguments are strings
    specifying the thermodynamic and kinetics model, respectively, and
    ``species`` and ``reactions`` keyword arguments are lists of `Species` and
    `Reaction` objects, respectively. Note that importing the reactions from a
    YAML input file requires a `Kinetics` object containing the species, as
    shown.

    Types of underlying models that form the composite `Solution` object are
    queried using the ``thermo_model``, ``kinetics_model`` and
    ``transport_model`` attributes; further, the ``composite`` attribute is a
    shorthand returning a tuple containing the types of the three constitutive
    models.

    For non-trivial uses cases of this functionality, see the examples
    `extract_submechanism.py <https://cantera.org/examples/python/kinetics/extract_submechanism.py.html>`_
    and `mechanism_reduction.py <https://cantera.org/examples/python/kinetics/mechanism_reduction.py.html>`_.

    In addition, `Solution` objects can be constructed by passing the text of
    the YAML phase definition in directly, using the ``yaml`` keyword
    argument::

        yaml_def = '''
        phases:
        - name: gas
          thermo: ideal-gas
          kinetics: gas
          elements: [O, H, Ar]
          species:
          - gri30.yaml/species: all
          reactions:
          - gri30.yaml/reactions: declared-species
          skip-undeclared-elements: true
          skip-undeclared-third-bodies: true
          state: {T: 300, P: 1 atm}
        '''
        gas = ct.Solution(yaml=yaml_def)
    """
    __slots__ = ()
    
    
    
def main():
    
    gas= Solution('gri30.yaml')
    
    return
    
if __name__ == '__main__': 
        main()    
        plt.show()
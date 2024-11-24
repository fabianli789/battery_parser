from typing import (
    TYPE_CHECKING,
)

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from structlog.stdlib import (
        BoundLogger,
    )

import numpy as np
from nomad.config import config
from nomad.datamodel.data import Schema
from nomad.metainfo import Quantity, SchemaPackage, Section, MSection, SubSection
from runschema.run import Run
from runschema.calculation import Calculation

configuration = config.get_plugin_entry_point(
    'battery_parser.schema_packages:schema_package_entry_point'
)

m_package = SchemaPackage()

class Dimensions(MSection):
    m_def = Section(validate=False)
    
    x = Quantity(type=float, description='x-dimension of 2d lattice in nm')
    y = Quantity(type=float, description='y-dimension of 2d lattice in nm')
    thickness = Quantity(type=float, description='final thickness of SEI in nm')

class Concentrations(MSection):
    m_def = Section(validate=False)
    
    name = Quantity(type=str, description='name of molecule')
    concentration = Quantity(type=np.float64, shape = ['*'], description='concentration of molecule at specific time, same length as "time"-array at run.calculation.concentration_time')


class Chem_Reactions(MSection):
    m_def  = Section(validate=False)
    
    name = Quantity(type=str, description = 'name of chem. reaction')
    barrier = Quantity(type=float, shape=[], description='energetic barrier in eV')
    escaped = Quantity(type=int, shape=[], description='number of escaped species')
    occurences = Quantity(type=np.float64, shape=[], description ='number of occurences of this chem. reaction')
    residence_time = Quantity(type=float, shape=[], description =  'time of each chem. reaction')
    

class BatteryCalculation(Calculation):
    m_def = Section(validate=False, extends_base_section=False)    
    
    dimensions = SubSection(sub_section=Dimensions.m_def, repeats=False)
    chem_reactions = SubSection(sub_section=Chem_Reactions, repeats=True)
    concentrations = SubSection(sub_section=Concentrations, repeats=True)
    volume_fraction = Quantity(type=float, description='volume if SEI got pressed together')
    porosity = Quantity(type=float, description='share of porous volume wrt to total SEI volume')
    concentration_time = Quantity(type=np.float64, shape=['*'], description='time evolution for the concentration of molecules, same length as "concentration"-array under run.calculation.concentrations')
    molecule_positions = Quantity(type=np.float64, shape=['*', 3], description='2D cartesian coordinates of molecules in nm.')
    molecule_species = Quantity(type=str, shape=['*'], description='Molecule species, same array length as molecule_positions.')


m_package.__init_metainfo__()

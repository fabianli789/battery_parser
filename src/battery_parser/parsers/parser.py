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

import yaml
import os
import re
import datetime
import numpy as np
import filecmp
from pathlib import Path

from nomad.config import config
from nomad.datamodel.metainfo.workflow import Workflow
from nomad.parsing.parser import MatchingParser
from runschema.run import Run, Program
from runschema.calculation import Calculation
from battery_parser.schema_packages.schema_package import BatteryCalculation, Chem_Reactions, Concentrations, Dimensions

configuration = config.get_plugin_entry_point(
    'battery_parser.parsers:parser_entry_point'
)



def DetailedBatteryParser(filepath, archive):
    run = Run()
    archive.run.append(run)
    run.program = Program(name="Meysam Battery Parser")
    
    calculation = BatteryCalculation()
    run.calculation.append(calculation)

    '''
    with open(str(filepath.parent) + r'/status_battery.csv') as status_file:
        time_run = archive.m_setdefault("run.time_run")
        time_run.cpu1_start = 0
        calc = archive.m_setdefault("run.calculation")
        _cpu = 0 
        _time = 0
        _step = 0
        for i, line in enumerate(status_file):
            line = line.strip("\n")
            parts = line.split(",")
            if parts[0] == None:
                continue
            if 'CPU' in line :
                try:
                    _cpu = float(parts[1])
                    time_run.cpu1_end = _cpu
                except:
                    _cpu = float('nan')
                    time_run.cpu1_end = _cpu
            if 'KMC time' in line:
                try:
                    _time = float(parts[1])
                    calc.time = _time
                except:
                    _time = float('nan')
                    calc.time = _time
            if 'steps' in line:
                try:
                    _step = int(float(parts[1]))
                    calc.step = _step
                except:
                    _step = -1
                    calc.step = _step
    '''
    with open(str(filepath.parent) + r'/SEI_properties_battery.csv') as prop_file:
        
        dimensions = Dimensions()
        calculation.dimensions = dimensions
        

        for i, line in enumerate(prop_file):
            line = line.strip("\n")
            parts = line.split(",")
            if parts[0] == None:
                continue
            if re.search(r'Thickness', parts[0]):
                dimensions.thickness = float(parts[1])    
            if re.search(r'Volume', parts[0]):
                calculation.volume_fraction = float(parts[1])
            if re.search(r'Poro', parts[0]):
                calculation.porosity = float(parts[1])

    with open(str(filepath.parent) + r"/occurrence_res_battery.csv") as occurrence_file:
        occurence_array = []
        residence_time_array = []
        for i, line in enumerate(occurrence_file):
            if re.search("Occurrences", line):
                continue
            part1, part2, part3 = line.split(",")
            try:
                occurence_array.append(int(part2))
            except:
                part2 = float('nan')
                occurence_array.append(part2)
            try:
                residence_time_array.append(float(part3))
            except:
                part3 = float('nan')
                residence_time_array.append(part3)
    with open(str(filepath.parent) + r'/concentration_battery.csv') as conc_file:
        first_line_parts = conc_file.readline().strip("\n").split(",")
        rows = 0
        for x, bla in enumerate(conc_file):
            rows = x
        conc_array = np.zeros((rows+1, len(first_line_parts)-2))
        time_array = []

        conc_file.seek(0) 
        for j, line in enumerate(conc_file):
    
            if re.search(r'time', line):
                continue
            
            parts = line.strip("\n").split(",")
            parts = [float(x) for x in parts]
            time_array.append(parts[-1])
            parts = parts[1:-1]            
            conc_array[j-1] = parts

        if rows > 0:
            calculation.concentration_time = time_array
        
            for i in range(len(first_line_parts)-2):
                
                concentrations = Concentrations()
                calculation.concentrations.append(concentrations)
                   
                concentrations.name = first_line_parts[i+1]
                concentrations.concentration = conc_array[:,i] 
        else:
            pass
    with open(str(filepath.parent) + r'/input_battery.yml') as file:
            j = 0
            for i, line in enumerate(file):
                parts  = line.split(": ")
                
                if "T" in parts[0]:
                    calculation.temperature = float(parts[1])
                if "xdim" in parts[0]:
                    dimensions.x = float(parts[1])
                if "ydim" in parts[0]:
                    dimensions.y = float(parts[1])
                if re.search(r'\-', parts[0]):
                    chem_reactions = Chem_Reactions()
                    calculation.chem_reactions.append(chem_reactions)                    
                    parts[0] = parts[0].lstrip('- ').rstrip(' ')
                    chem_reactions.name = parts[0]
                    if re.search(r'escape', parts[0]):
                            
                        escaped = Escaped(Path(filepath).parent, chem_reactions)
                        chem_reactions.escaped = escaped[j]
                        j += 1
                    chem_reactions.barrier = float(parts[1])
                    
                    chem_reactions.occurences = occurence_array[i-5]
                    chem_reactions.residence_time = residence_time_array[i-5]
    with open(str(filepath.parent) + r'/last_step_battery.csv') as last_step_file:
        species_array = []
        for j, x in enumerate(last_step_file):
            pass
        coordinates = np.zeros((j, 3))
        coord_x = []
        coord_y = []
        last_step_file.seek(0)    
        for i, line in enumerate(last_step_file):
            
            parts = line.strip("\n").split(",")
            
            if re.search(r'species', line):
                continue
            try:
                coord_x.append(float(parts[0].strip('"').strip("[")))
            except:
                coord_x.append(float('nan'))
            try:
                coord_y.append(float(parts[1].strip('"').strip("]")))
            except:
                coord_y.append(float('nan'))

            species_array.append(parts[2])
        coordinates[:, 0] = coord_x
        coordinates[:, 1] = coord_y
        calculation.molecule_positions = coordinates
        
        calculation.molecule_species = species_array
def Escaped(parent, chem_reactions):
    escaped_file =str(parent)+r"/escaped_battery.csv"
    with open(escaped_file) as file:
        escaped_array = []
        for i, line in enumerate(file):
            parts =  line.strip("\n").split(",")
            if len(parts) ==  0 or len(parts) == 1:
                continue
            
            escaped_array.append(int(parts[1])) 
    return escaped_array        



class NewParser(MatchingParser):
    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] = None,
    ) -> None:
        logger.info('NewParser.parse', parameter=configuration.parameter)

        archive.workflow2 = Workflow(name='test')
        
        mainfile = Path(mainfile)
        DetailedBatteryParser(mainfile, archive)
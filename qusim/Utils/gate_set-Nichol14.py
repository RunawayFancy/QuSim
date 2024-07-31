from qusim.PulseGen.pulse_config import PulseConfig
from copy import deepcopy
import json
import os
from collections import defaultdict as ddict
from typing import Literal, Optional
from tkinter import Tk, filedialog
import numpy as np

__SINGLE_Q_GATE_LABEL__ = [
    "I", "X", "Y", "Z", 'H'
]
__TWO_Q_GATE_LABEL__ = [
    "CZ", "CNOT", "SWAP" 
]

__DEFAULT_PATH__ = "gate_set/"

# This might do the same work comparing to smdata...

class SystemInfo:
    """
    t_moment: float
        Referring to https://quantumai.google/cirq/build/circuits
    
    """
    def __init__(self, num_qubit: int , t_moment: float) -> None:
        self.num_qubit = num_qubit
        self.t_moment = t_moment
    
    @property
    def sysinfo2dict(self) -> dict:
        system_info_dict = dict()
        system_info_dict["num_qubit"] = self.num_qubit
        system_info_dict["t_moment"] = self.t_moment

        return system_info_dict


class GateSet:
    def __init__(
            self,
            path: str,
            case: Literal["load", "new"],
            filename: Optional[str] = None,
            system_info: Optional[SystemInfo] = None
        ) -> None:

        if path.endswith("/"):
            self.path = f"{path}{__DEFAULT_PATH__}"
        else:
            self.path = f"{path}/{__DEFAULT_PATH__}"
        self.init_directory
        self.case = case
        if self.case == "new":
            if not system_info:
                raise(ModuleNotFoundError("Missing system_info for create a new GateSet."))
            if not filename:
                raise(ModuleNotFoundError("Missing filename for create a new GateSet."))
            self.system_info = system_info
            self.filename = filename
            self.gate_param = ddict(list)
            self.gate_info = dict()
            self.gate_set_json = {"system_info": self.system_info.sysinfo2dict, "gate_info": self.gate_info, "gate_param": self.gate_param}
        else:
            if not filename:
                self.filename = self.sfile
            else:
                self.filename = filename
            self.gate_set_json = json.load(open(self.path+self.filename, 'rb'))
            self.system_info = SystemInfo(self.gate_set_json["system_info"]["num_qubit"], self.gate_set_json["system_info"]["t_moment"])
            self.gate_info = self.gate_set_json["gate_info"]
            self.gate_param = self.gate_set_json["gate_param"]
        pass


    def get_num_of_moment(self, gate_label: str):
        try:
            gate_param_lst = self.gate_param.get(gate_label,0)
        except:
            raise(ValueError(f"The input {gate_label} should be a string."))
        if gate_param_lst ==0:
            raise(ValueError(f"Not such gate_label {gate_label} is found in gate_param dict."))
        t_delay_lst = [pulse_config_dict["t_delay"] for pulse_config_dict in gate_param_lst]
        num_of_moment = len(np.unique(t_delay_lst))

        return num_of_moment


    def gate_param_validity(self, pulse_param_dict):
        assert pulse_param_dict["t_width"] + pulse_param_dict["t_plateau"] <= self.system_info.t_moment
        pass


    def add(self, gate_label: str, pseq: list[PulseConfig]):
        """
        gate_label: str
            The gate_label should take the form "GATE_Q1&Q2&Q3". 
            E.g., a Pi rotation along X axis on qubit 0 (qindex 0) should be: "X90_0";
            A Pi/2 rotation with phase 34 degree on qubit 2 (qindex 2) should be "R90P34_0"
            A identity gate on qubit 1 (qindex 1) should be "I_1"
            A Z rotation gate with angle 34.1 degree on qubit 5 (qindex 6) should be "Z34_5"
            A CZ gate between qubit 1 (qindex 1) and qubit 2 (qindex 2) is "CZ_1&2"
        """
        switch = 1
        if gate_label in self.gate_param.keys():
            response = input(f"Gate parameter already in gate_param dict, still overwrite it? [y/n]: ")
            if response.lower() == 'y':
                pass
            elif response.lower() == 'n':
                switch = 0
                print("No action taken, think twice next time!")
            else:
                switch = 0
                print("Invalid input. Please enter 'y' or 'n'. ")
        if switch > 0:
            gate_param_lst_buffer = []
            for pulse in pseq:
                pulse_param_dict = pulse.export2dict()
                self.gate_param_validity(pulse_param_dict)
                gate_param_lst_buffer.append(pulse_param_dict)
            self.gate_param[gate_label] = gate_param_lst_buffer

            t_delay_lst = [[_i, pulse_config_dict["t_delay"]] for _i,pulse_config_dict in enumerate(self.gate_param[gate_label])]
            sorted_t_delay_lst = sorted(t_delay_lst, key=lambda x: x[1])

            current_t_delay = np.min(np.array(sorted_t_delay_lst)[:,1])
            _j=0
            for kv in sorted_t_delay_lst:
                if kv[1] > current_t_delay:
                    _j += 1
                self.gate_param[gate_label][kv[0]]["t_delay"]=self.system_info.t_moment * _j

            num_of_moment = len(np.unique(t_delay_lst))
            try:
                self.gate_info[gate_label]["num_of_moment"] = num_of_moment
            except:
                self.gate_info[gate_label] = ddict()
                self.gate_info[gate_label]["num_of_moment"] = num_of_moment
        pass


    def delete(self, gate_label: str):
        response = input(f"Deleting gate parameter {gate_label}, are you sure? [y/n]: ")
        if response.lower() == 'y':
            del self.gate_param[gate_label]
        elif response.lower() == 'n':
            print("No action taken, think twice next time!")
        else:
            print("Invalid input. Please enter 'y' or 'n'. ")
        pass
    

    @property
    def sleep(self):
        print("Sleeping...")
        json.dump(self.gate_set_json, open(f"{self.path}{self.filename}.json", "w"), sort_keys=True, separators=(',', ': '))
        print(f"File save as {self.path}{self.filename}")
        pass


    @property
    def init_directory(self):
        if not os.path.exists(self.path):
            try:
                os.makedirs(self.path)
                print(f'Directory {self.path} created successfully.')
            except Exception as e:
                print(f"An error occurred while creating directory {self.path}: {e}")
        else:
            print(f'Directory {self.path} exists.')


    @property
    def reload(self):
        response = input(f"Reloading gate parameter, are you sure? [y/n]: ")
        if response.lower() == 'y':
            self.gate_set_json = json.load(open(f"{self.path}{self.filename}.json", 'rb'))
            system_info_dict = self.gate_set_json["system_info"]
            self.system_info = SystemInfo(system_info_dict["num_qubit"], system_info_dict["t_moment"])
            self.gate_param = self.gate_set_json["gate_param"]
        elif response.lower() == 'n':
            print("No action taken, think twice next time!")
        else:
            print("Invalid input. Please enter 'y' or 'n'. ")
        pass


    @property
    def sfile(self):
        root = Tk()
        root.withdraw()
        root.call('wm', 'attributes', '.', '-topmost', True)
        filepaths = filedialog.askopenfilenames(initialdir=self.path)
        root.destroy()
        filenames = [os.path.basename(filepath) for filepath in filepaths]
        assert len(filenames)==1 # Only one file can be selected
        print(f"Selected files: {filenames}")

        return filenames  # Return the list of filenames for further use


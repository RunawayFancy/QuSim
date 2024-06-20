# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""
import os
import pickle
from tkinter import Tk, filedialog

class qsave:

    def __init__(self, path):
        self.path = path
        self.init_directory()

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
    def init(self):
        try:
            counts = pickle.load(open(self.path+'counts.pkl', 'rb'))
            response = input("Current counts.pickle exist, do you want to reinitialize the counts.pkl file? [y/n]: ")
            if response.lower() == 'y':
                switch = 1
            elif response.lower() == 'n': 
                switch = 0
                print("No action taken.")
            else: print("Invalid input. Please enter 'y' or 'n'.")
        except FileNotFoundError:
            switch = 1
        if switch == 1:
            counts = [0, '']
            pickle.dump(counts, open(self.path+'counts.pkl', 'wb'))
            print(f'Successfully create counts.pickle file under directory {self.path}')

    @property
    def check(self):
        try:
            counts = pickle.load(open(self.path+'counts.pkl', 'rb'))
            print(f"Current counts.pickle exist, with value: {counts}")
        except FileNotFoundError:
            print(f"Current counts.pickle missing, please initial a new one using `qsave.init`")

    @property
    def delete(self):
        response = input("Do you want to delete the counts.pkl file? [y/n]: ")
        if response.lower() == 'y':
            try:
                os.remove(self.path + 'counts.pkl')
                print(f'Successfully deleted counts.pkl file under directory {self.path}')
            except FileNotFoundError:
                print(f"The file {self.path}counts.pkl does not exist.")
            except PermissionError:
                print(f"Permission denied: unable to delete {self.path}counts.pkl.")
            except Exception as e:
                print(f"An error occurred: {e}")
        elif response.lower() == 'n':
            print("No action taken.")
        else:
            print("Invalid input. Please enter 'y' or 'n'.")

    @property
    def undo(self):
        count, filename = pickle.load(open(self.path+"counts.pkl", "rb"))
        # print("Undoing the process...")
        fpath = self.path+filename
        response = input(f"Do you want to undo the {filename}? [y/n]: ")
        # Check the user's input and respond
        if response.lower() == 'y':
            try:
                os.remove(fpath)
                print(f"Successfully undo {filename}")
            # Add code here to undo the process
            except FileNotFoundError:
                print(f"The file {fpath} does not exist.")
            except PermissionError:
                print(f"Permission denied: unable to delete {fpath}.")
            except Exception as e:
                print(f"An error occurred: {e}")
            if count > 0:
                count -= 1

        elif response.lower() == 'n':
            print("No action taken.")
        else:
            print("Invalid input. Please enter 'y' or 'n'.")
        pickle.dump([count, filename], open("counts.pkl", "wb"))

    def save(self, filename: str ,data=None, scan=None):
        count, _ = pickle.load(open(self.path + "counts.pkl", "rb"))
        filename = filename + f'_{count}.pkl'
        # scan save, need a formate scan dict.
        if scan is not None:
            exit(0)
        # custom data save
        if data is not None:
            pickle.dump(data, open(self.path+filename, "wb"))
            pickle.dump([count+1, filename], open(self.path+"counts.pkl", "wb"))
            print("Save as file "+filename)

    @property
    def sfile(self):
        root = Tk()
        root.withdraw()
        root.call('wm', 'attributes', '.', '-topmost', True)
        filepaths = filedialog.askopenfilenames(initialdir=self.path)
        root.destroy()
        filenames = [os.path.basename(filepath) for filepath in filepaths]
        print("Selected files:", filenames)
        return filenames  # Return the list of filenames for further use
    
    @property
    def load(self):
        filenames = self.sfile
        data = []
        for fn in filenames:
            data.append(pickle.load(open(self.path + fn, 'rb')))
        return data

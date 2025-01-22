import os
from pathlib import Path
import pickle
import lzma


import sys
sys.stdout.flush()
try:
    shell = get_ipython().__class__.__name__
    if shell == 'ZMQInteractiveShell': # script being run in Jupyter notebook
        from tqdm.notebook import tqdm
    elif shell == 'TerminalInteractiveShell': #script being run in iPython terminal
        from tqdm import tqdm
except NameError:
    if sys.stderr.isatty():
        from tqdm import tqdm
    else:
        from tqdm import tqdm # Probably runing on standard python terminal. If does not work => should be replaced by tqdm(x) = identity(x)

parent_dir =  'AddYourDefaultFolderHERE'

def save_preprocess(filename, saving_type, parent_dir):
    """Creates parent dir if does not exist and ensures filename has the saving_type extension."""
    Path(parent_dir).mkdir(parents=True, exist_ok=True)
    filename = filename if filename.endswith('.{}'.format(saving_type)) else '{}.{}'.format(filename, saving_type)
    return filename
        
def create_saving_func(file_opener, file_writer, saving_type, **kwargs):
    def saving_func(file, filename, parent_dir=parent_dir):
        filename = save_preprocess(filename, saving_type, parent_dir)
        with file_opener(os.path.join(parent_dir, filename), 'w') as f:
            file_writer(file, f, **kwargs)
        return
    return saving_func

save_lzma = create_saving_func(lzma.LZMAFile, pickle.dump, 'lzma')
        

def create_loading_func(file_opener, file_loader, extra_processing=None, apply_defaults=None):
    """
    - extra_processing = List of process that could be applied to the file.
    - apply_defaults = dict:{str: bool}. The key is set as a function variable, the value indicates whether to apply the process named by the key.
    """
    if extra_processing is None:
        def loading_func(path):
            with file_opener(path, 'rb') as f:
                return file_loader(f)
    else: 
        return None
        def loading_func(path):
            args = {key:val for key,val in locals().items() if key not in  ('path', *create_loading_func.__code__.co_varnames)}
            with file_opener(path, 'rb') as f:
                data = file_loader(f)
            
            for condition, process in zip(args.values(), extra_processing):
                if condition:
                    data = process(data)
            return data
        
        # I know this sould not be done this way, but I wanted to check it can
        code_data = loading_func.__code__
        num_vars = code_data.co_argcount
        num_vars_new = len(apply_defaults)
        new_code = code_data.replace(co_varnames=(*code_data.co_varnames[:num_vars], *apply_defaults.keys(), *code_data.co_varnames[num_vars:]),
                                     co_argcount=num_vars + num_vars_new, 
                                     co_nlocals=code_data.co_nlocals + num_vars_new)
        loading_func.__code__ = new_code
        loading_func.__defaults__ = tuple(apply_defaults.values())
    return loading_func

load_lzma = create_loading_func(lzma.LZMAFile, pickle.load)
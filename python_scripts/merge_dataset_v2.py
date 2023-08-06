import os
import numpy as np
from scipy.io import loadmat, savemat
from datetime import datetime


def file2List(directory, ending):
    merged_data = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(ending):
                filepath = os.path.join(root, file)            
                with open(filepath, 'r') as file:
                    data = file.read().replace('\n', ',')
                merged_data = merged_data + data.split(',')

    return merged_data



def merge_mat_files(directory, ending, merged_data=[], nmax=100):
    new_merged_data = []
    result = {}
    # Walk through the directory tree
    i = 0
    first_entry = True
    for root, _, files in os.walk(directory):
        for file in files:
            if ending in file: 
                filepath = os.path.join(root, file)
                # Get the keys
                if first_entry:
                    bound = loadmat(filepath)
                    for keys, _ in bound.items():
                        result[keys] = []
                    first_entry = False
                # Extract the file name without extension
                dirname = root[len(directory)+1::]
                # print(dirname)
                # Load the data from the MATLAB file and Merge the data into the dictionary if not yet merged
                if dirname not in merged_data:
                    i = i+1 
                    new_merged_data.append(dirname)
                    bound  = loadmat(filepath)
                    print("Saving folder %s"%(i))   
                    for keys, _ in bound.items():
                        result[keys].append(bound[keys])
                else:
                    continue     
                    
    for keys, _ in bound.items():
        if (keys in list_special):
            # print(keys)
            continue
        if result[keys]:    
            # print(keys)
            result[keys] = np.concatenate(result[keys], axis=0)
            #print(type(result[keys]))
            batch_size =   result[keys].shape[0]  
            # print(result[keys].shape)
            # print(batch_size)

    if new_merged_data:
        with open('%s/logs_mergedFiles_%s.%s'%(directory_path, TIMESTAMP, ending[:]), 'w') as f:
            for line in new_merged_data:
                f.write("%s\n" % line)
        
        for keyss in list_special: # These keys have values which are same for all batches, so do not need to be repeated
            if keyss in result:
                result[keyss]= result[keyss][0]

        start = 0
        stop = 0
        RESULT = {}
        i = 0
        while start < batch_size:
            stop = min(start + nmax, batch_size)  

            for keys, values in result.items():
                if (keys in list_special): 
                    #print(keys)
                    RESULT[keys] = values
                    #continue
                else:
                    RESULT[keys] = values[start:stop, ...]

            output_folder = '%s/%s_samples__max_Inclusions_3__%s__%s'%(directory_path, stop-start, TIMESTAMP, i)
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)

            savemat('%s/%s'%(output_folder, savename[ending]), RESULT)  
            start = start + nmax 
            i = i+1
            print("\nCopy ended at index %s"%(stop)) 
            #print(batch_size)

    return merged_data, result


if __name__ == "__main__":
    # Replace 'your_directory_path' with the path to your main directory
    TIMESTAMP = datetime.utcnow().strftime('%Y%m%d-%H%M%S-%f')
    savename = {}
    savename['bound' ] = 'dataset_bound.mat'
    savename['domain'] = 'dataset_domain.mat'
    savename['mesh'  ] = 'mesh.mat'
    list_special = ["angl_circum", "radius", "theta", "x1", "x2", "__header__", "__version__", "__globals__"]

    directory_path = '/pvfs2/Derick/EIT/Mine/data'

    for ending in ['bound']:#['domain', 'bound', ]:
        #ending = "bound" #"dataset_bound.mat"
        nmax = 10000

        merged_data = file2List(directory_path, '.txt')

        merged_data, result = merge_mat_files(directory_path, ending, merged_data, nmax)

    #print()





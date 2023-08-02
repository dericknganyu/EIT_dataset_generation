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
    # Dictionary to store the merged data
    
    new_merged_data = []
    result = {}

    # Walk through the directory tree
    i = 0
    first_entry = True
    
    for root, _, files in os.walk(directory):
        for file in files:
            if ending in file: #file.endswith(ending):
                filepath = os.path.join(root, file)

                if first_entry:
                    dataset = loadmat(filepath)
                    if 'outputVoltage' in dataset.keys():
                        dataset.pop('outputVoltage')
                    for keys, _ in dataset.items():
                        result[keys] = []
                    first_entry = False

                # Extract the file name without extension
                dirname = root[len(directory)+1::]
                # print(dirname)

                # Load the data from the MATLAB file
            
                # Merge the data into the dictionary
                if dirname not in merged_data:
                    #print('Yeah')
                    i = i+1 
                    new_merged_data.append(dirname)
                    dataset  = loadmat(filepath)
                    if 'outputVoltage' in dataset.keys():
                        dataset.pop('outputVoltage')
                    print("Saving folder %s : %s"%(dirname, i))   
                    for keys, _ in dataset.items():
                        result[keys].append(dataset[keys])
                else:
                    continue



    for keys, _ in dataset.items():
        if (keys in list_special):
            #print(keys)
            continue
        if result[keys]:    
            print(keys)
            result[keys] = np.concatenate(result[keys], axis=0)
            #print(type(result[keys]))
            batch_size =   result[keys].shape[0]  
            print(result[keys].shape)
            #print(batch_size)

    if new_merged_data:
        with open('%s/logs_mergedFiles_%s.%s'%(directory_path, TIMESTAMP, ending[:]), 'w') as f:
            for line in new_merged_data:
                f.write("%s\n" % line)
        
        for keyss in list_special:
            if keyss in result:
                result[keyss]= result[keyss][0]

        # for keys, _ in result.items():
        #     if keys == "__header__" or keys == "__version__" or keys == "__globals__":
        #         #print(keys)
        #         continue
        #     if result[keys].any():    
        #         print(result[keys].shape)


        start = 0
        stop = 0
        RESULT = {}
        i = 0
        while start < batch_size:
            stop = min(start + nmax, batch_size)  

            for keys, values in result.items():
                if (keys in list_special): #keys == "__header__" or keys == "__version__" or keys == "__globals__" or keys in list_special:
                    #print(keys)
                    RESULT[keys] = values
                    #continue
                else:
                    RESULT[keys] = values[start:stop, ...]

            output_folder = '%s/%s_samples__max_Inclusions_3__%s__%s'%(directory_path, stop-start, TIMESTAMP, i)
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
            #RESULT.pop('outputVoltage')
            savemat('%s/%s'%(output_folder, savename[ending]), RESULT)  
            start = start + nmax 
            i = i+1
            print("\nCopy ended at index %s"%(stop)) 
            print(RESULT.keys())
            print()


    return merged_data, result


if __name__ == "__main__":
    # Replace 'your_directory_path' with the path to your main directory
    TIMESTAMP = datetime.utcnow().strftime('%Y%m%d-%H%M%S-%f')
    savename = {}
    savename['bound' ] = 'dataset_bound.mat'
    savename['domain'] = 'dataset_domain.mat'
    savename['mesh'  ] = 'mesh.mat'
    list_special = ["angl_circum", "radius", "theta", "x1", "x2", "__header__", "__version__", "__globals__"]

    for j in range(2,14):
        directory_path = '/pvfs2/Derick/EIT/Mine/data_%s'%(j)
        print("\nWorking on %s"%(directory_path))
        for ending in ['domain', 'bound']:
            print("\nWorking on %s"%(savename[ending]))
            #ending = "bound" #"dataset_bound.mat"
            nmax = 5000

            merged_data = file2List(directory_path, '.txt')

            merged_data, result = merge_mat_files(directory_path, ending, merged_data, nmax)

    #print()





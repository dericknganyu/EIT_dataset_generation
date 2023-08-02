from pydrive2.auth import GoogleAuth
from pydrive2.drive import GoogleDrive
import os

import time

def upload_folder_to_drive(folder_path, parent_id=None):
    gauth = GoogleAuth()
    gauth.LocalWebserverAuth()

    drive = GoogleDrive(gauth)

    folder_name = os.path.basename(folder_path)

    folder_metadata = {'title': folder_name, 'mimeType': 'application/vnd.google-apps.folder'}
    if parent_id:
        folder_metadata['parents'] = [{'id': parent_id}]

    folder = drive.CreateFile(folder_metadata)
    folder.Upload()
    i = 0
    for item in os.listdir(folder_path):
        item_path = os.path.join(folder_path, item)

        if os.path.isdir(item_path):
            print('\nUploading directory >> %s folder %s'%(item_path, i+1))
            upload_folder_to_drive(item_path, parent_id=folder['id'])
        else:
            #if item_path.endswith('.png') and not('400' in item_path):
            print(  'Uploading file      >> %s'%(item_path))
            file = drive.CreateFile({'title': item, 'parents': [{'id': folder['id']}]})
            file.SetContentFile(item_path)
            file.Upload()
            
        i = i+1

if __name__ == "__main__":
    # Replace 'YOUR_FOLDER_PATH' with the path to the directory you want to upload.
    folder_path = '/pvfs2/Derick/EIT/Mine/dataset' #/pvfs2/Derick/EIT/Mine/data'#'/localdata/Derick/EIT/Mine/data'
    #ID = "1qwcoVE5VTQ0AQzO3LTtMzRBblGy59KwL" #moriarty
    #ID = "1P0MMqQf7Juf5-7-TLOAe5yVyyxH2_UnI" #data jarvis
    #ID = "19nDaXvJlXzOm-edUAn5FX2KCf-Muw6C8" #data ztm-cluster
    ID = "1zCqWx26ZNWiYRWjI-Bd7JvVDo8zmr6GP" #dataset

    t0 = time.time()

    upload_folder_to_drive(folder_path, ID)

    print("_______________________________________________________")
    print("Process ends after %.4f minutes"%((time.time()-t0)/60))
    print("_______________________________________________________")

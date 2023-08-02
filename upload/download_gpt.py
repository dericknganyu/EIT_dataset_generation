from pydrive2.auth import GoogleAuth
from pydrive2.drive import GoogleDrive

# Function to download a file or folder from Google Drive recursively
def download_folder_from_drive(folder_id, output_folder):
    # Authenticate and create GoogleDrive instance
    gauth = GoogleAuth()
    gauth.LocalWebserverAuth()
    drive = GoogleDrive(gauth)

    # Retrieve the folder information
    folder = drive.CreateFile({'id': folder_id})
    folder.FetchMetadata(fields='title')

    print("Downloading folder: %s"%(folder['title']))

    # Create the output folder if it doesn't exist
    import os
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Retrieve the list of files and folders inside the main folder
    file_list = drive.ListFile({'q': "'{}' in parents and trashed=false".format(folder_id)}).GetList()

    # Download each file/folder recursively
    for file_or_folder in file_list:
        if file_or_folder['mimeType'] == 'application/vnd.google-apps.folder':
            # If it's a sub-folder, recursively download it
            subfolder_path = os.path.join(output_folder, file_or_folder['title'])
            download_folder_from_drive(file_or_folder['id'], subfolder_path)
        else:
            # If it's a file, download it
            file_path = os.path.join(output_folder, file_or_folder['title'])
            file_or_folder.GetContentFile(file_path)

if __name__ == "__main__":
    # Replace 'YOUR_FOLDER_ID' with the ID of the folder you want to download
    folder_id = '12ZxGXaBx4eZBNqlPY_LJLFH7z6-al-Ux' #moriarty
    #folder_id = '1P0MMqQf7Juf5-7-TLOAe5yVyyxH2_UnI' #jarvis
    
    # Replace 'YOUR_OUTPUT_FOLDER_PATH' with the local path where you want to save the downloaded folder
    output_folder = '/pvfs2/Derick/EIT/Mine/jarvis' #moriarty

    download_folder_from_drive(folder_id, output_folder)

    # 1P0MMqQf7Juf5-7-TLOAe5yVyyxH2_UnI
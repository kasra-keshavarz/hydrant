

def merge_multiple_keys(land_cover = 'CEC',
                        keys = ['1','2']):
    
    from .MESH_default_dict import LULC
    
    merged_data = ""
    for key in keys:
        data = LULC.get(key)
        if data:
            merged_data += f"--- Key {key} ---\n{data}\n\n"

    # if merged_data:
    #     with open("merged_data.txt", "w") as file:
    #         file.write(merged_data)
    #     print(f"Merged data for keys {', '.join(keys)} has been written to merged_data.txt")
    # else:
    #     print("No data found for the given keys")

# Example usage:
merge_multiple_keys(['1', '2', '3'])  # Pass the list of keys you want to merge
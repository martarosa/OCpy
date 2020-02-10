def check_namelist_key_and_print(key, output_string):
    if check_namelist_key_exist(key):
        print(output_string)


def check_namelist_key_exist_and_value(key, value):
    if key in locals() or key in globals():
        if key == value:
            return True


def check_namelist_key_exist_and_list_value(key, value):
    if key in locals() or key in globals():
        for i in range(len(value)):
            if key == value[i]:
                return True

def check_namelist_key_exist(key):
    if key in locals() or key in globals():
        return True


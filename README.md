```shell
# show pip version (and check it is installed)
py -m pip --version

# show installed libraries
py -m pip freeze 

# write installed libraries to requirements.txt file
py -m pip freeze >requirements.txt

# reinstall libraries (e.g. on another machine)
py -m pip install -r requirements.txt
```
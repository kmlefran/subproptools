Tips
====

- Units are in atomic units throughout

- When calling functions, refer to atoms by the number in the .log files, not python Indices
    - That is, start at 1, not 0

- Some functions require their input in the form of a list of strings of the lines of a .sum file. This input can be generated using code similar to the following:

.. code-block:: python

    with open('path/to/file.sum','r',encoding='utf-8') as file:
        data = file.readlines()

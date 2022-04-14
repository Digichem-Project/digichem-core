Installation
============

To install Silico, connect to Kennedy using PuTTY and log-in to your account as normal.
Once logged in, use the following command:

.. code-block:: console
	
	$ /gpfs1/apps/EZC-tools/install-silico
	
.. note::
	Installation only ever needs to be performed once (for each user). Updates to Silico will be become available automatically without the need to reinstall.
	
Installation should be near instantaneous. Once complete, a message will be printed:

.. code-block:: console

	$ /gpfs1/apps/EZC-tools/install-silico
	Installed successfully. Please log out and in again to complete
	
Use the ‘exit’ command to log out:

.. code-block:: console

	$ exit
	logout
	
Once logged back in again, Silico will now be available for use.
The installation can be tested by using the ``silico`` command with the ‘-v’ (version) option.
If successful, the version will be printed:

.. code-block:: console

	$ silico -v
	1.0.0

If the installation was not successfully, the following error message will be printed:

.. code-block:: console

	$ silico -v
	-bash: silico: command not found
	
If unsuccessful, seek help from another group member.
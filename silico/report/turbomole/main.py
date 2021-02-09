from silico.report.base.pdf import PDF_report


class Turbomole_report(PDF_report):
    """
    A specialised report object for processing Turbomole results.
    """
    
    def setup_cubes(self, output_dir, output_name):
        """
        Setup the cube files which will be used to render images.
        
        :param output_dir: A pathlib Path object to the directory within which our files should be created.
        :param output_name: A string that will be used as the start of the file name of the files we create.
        """
        # First, get our fchk file (from which cubes are made in Gaussian.
        self.fchk_file = Chk_to_fchk(
            Path(output_dir, output_name + ".fchk"),
            chk_file = self.chk_file_path,
            fchk_file = self.fchk_file_path,
        )
        
        
        
        ################
        # Spin density #
        ################
        self.cubes['spin_density'] = Fchk_to_spin_cube.from_options(
            Path(output_dir, "Spin Density", output_name + ".spin.cube"),
            fchk_file = self.fchk_file,
            spin_density = "SCF",
            options = self.options
        )
        
        
        ############
        # Orbitals #
        ############
        # We need to set images for both alpha and beta orbitals (if we have them).
        for orbital_list in (self.result.molecular_orbitals, self.result.beta_orbitals):
            for orbital in orbital_list:
                # First, decide what type of orbital we need.
                if orbital.spin_type == "alpha":
                    cubegen_type = "AMO"
                elif orbital.spin_type == "beta":
                    cubegen_type = "BMO"
                else:
                    cubegen_type = "MO"
                
                # Save cube.
                self.cubes[orbital.label] = Fchk_to_cube.from_options(
                    Path(output_dir, orbital.label, output_name + ".{}.cube".format(orbital.label)),
                    fchk_file = self.fchk_file,
                    cubegen_type = cubegen_type,
                    orbital = orbital.level,
                    options = self.options)
        
        
        #############
        # Structure #
        #############
        # If we have an orbital cube, we can just reuse this for our structure.
        if "HOMO" in self.cubes:
            self.cubes['structure'] = self.cubes['HOMO']
        elif "HOMO (alpha)" in self.cubes:
            self.cubes['structure'] = self.cubes['HOMO (alpha)']
        else:
            # No MO cubes available, create one for structure.
            # We'll just use the HOMO to get our cube, as it almost certainly should exist.
            self.cubes['structure'] = Fchk_to_cube.from_options(
                Path(output_dir, "Structure", output_name + ".structure.cube"),
                fchk_file = self.fchk_file,
                cubegen_type = "MO",
                orbital = "HOMO",
                options = self.options)
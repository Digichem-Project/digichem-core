Changes from the previous version
=================================

Version V1.0.0 is the first release of the V1.x branch of Silico. This release brings a swathe of new features and improvements, including a general overhaul of certain internals to improve the overall stability and performance of the program, as well as the first stable (ie, not liable to sudden change) interface. Notable changes from the previous, V0.x, build are as follows:

 * A fully-functioning interactive (urwid) interface to all of the major subprograms, and a reworked and improved interface to the submit subprogram.
 * Various changes to the console interfaces to all subprograms, particularly towards improving compatibility between the subprograms.
 * Introduced the ability to submit calculations from method files (previously only methods from the library could be used).
 * Added support for rendering natural-transition orbital plots (Gaussian).
 * Added support for rendering difference density plots (Turbomole).
 * Added a new report style (the journal style) that emulates a scientific journal article.
 * Entirely reworked how method definitions are handled in the internal library, improving performance and allowing a significantly greater number of definitions.
 * Various bug-fixes and minor improvements.

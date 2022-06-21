#import bibtexparser
from bibtexparser.bparser import BibTexParser
import bibtexparser.customization
from datetime import datetime

# Although we are using a bibtex parser for reading in our references, it's not that great and we're still forced to do a lot of work ourselves :(.
# TODO: Find a better bibtex library that can do more of what we want. 

class RSC_format():
    """
    Top-level class for printing references in RSC (Royal Society of Chemistry) format. 
    """
    def __init__(
        self,
        record
    ):
        """
        Constructor for Reference objects.
        """
        self.record = record
        
    @classmethod
    def bibtex_customizations(self, record):
        """
        This function is used by the bibtexparser parser to convert basic string types into more useful types.
        """
        # Splits the 'authors' field into a list of strings. The docs suggest this string is “Name, Surname” but I suspect it is in fact "Surname, Name".
        record = bibtexparser.customization.author(record)
        record = bibtexparser.customization.page_double_hyphen(record)
        return record
    
    @classmethod
    def from_record(self, record):
        """
        Get an appropriate RSC_format object depending on the type of a given record.
        """
        if record['ENTRYTYPE'] == "article":
            return RSC_format_journal(record)
        elif record['ENTRYTYPE'] == 'online':
            return RSC_format_online(record)
        elif record['ENTRYTYPE'] == 'software':
            return RSC_format_software(record)
        elif "thesis" in record['ENTRYTYPE']:
            return RSC_format_thesis(record)
        else:
            return self(record)
    
    @classmethod
    def dict_from_file(self, file_object):
        """
        Get a dictionary of RSC_format references (where each key is the unique name of the reference) from an open() .bib file.
        
        :param file_object: An open() file-like bibtex object to read from.
        :return: A dictionar of RSC_format objects.
        """
        # First load our list of references.
        parser = BibTexParser()
        parser.customization = self.bibtex_customizations
        parser.ignore_nonstandard_types = False
        bib_db = bibtexparser.load(file_object, parser=parser)
        
        # Now make a dictionary of RSC_format objects from the dict of references.
        return {key:self.from_record(reference) for key,reference in bib_db.entries_dict.items()}
    
    @classmethod
    def format_name(self, name):
        """
        Format a name into RSC style.
        
        :param name: The name as a string (as formatted by bibtexparser.customization.author.
        :returns: The name in RSC style (as a string).
        """
        # This format (von and Jr) might not be quite right...
        # First, split the name.
        name = bibtexparser.customization.splitname(name)
        # Abbreviate all first names.
        formatted_name = " ".join(["{}.".format(first_name[0].upper()) for first_name in name['first']])
        # Add the 'von' part if we have it
        if len(name['von']) > 0:
            formatted_name += " " + " ".join(name['von'])
        # Add surnames
        if len(name['last']) > 0:
            formatted_name += " " + " ".join(name['last'])
        # Add Jr.
        if len(name['jr']) > 0:
            formatted_name += ", " + " ".join(name['jr'])
        # All done.
        return formatted_name
    
    @property
    def authors(self):
        """
        Get the names of all the authors of this reference as a string.
        """
        unformatted_authors = self.record.get('author')
        if unformatted_authors is None:
            return None
        
        # First get a list of each author's name formatted correctly.
        names = [self.format_name(name) for name in unformatted_authors]
        
        # Now stitch them together.
        if len(names) == 1:
            return names[0]
        elif len(names) > 1:
            return " and ".join((", ".join(names[:-1]), names[-1]))
        else:
            return None
        
    @property
    def year(self):
        """
        The year of publication of this reference as a string.
        """
        return self.record.get('year')
    
    @property
    def url(self):
        """
        The URL of this reference as a string.
        """
        return self.record.get('url')
    
    @property
    def sections(self):
        """
        Get a list of the sections that make up a reference of this style.
        """
        return [self.authors, self.year]
    
    def __str__(self):
        """
        Stringify this reference.
        """
        # First get a list of sections that aren't empty.
        sections = [ref_part for ref_part in self.sections if ref_part is not None]
        # Return the none empty sections tuck together.
        return ", ".join(sections)
    

class RSC_format_journal(RSC_format):
    """
    Class for printing journal references in RSC (Royal Society of Chemistry) format.
    """    
        
    @property
    def journal(self):
        """
        The journal of this reference as a string.
        """
        # TODO: This should use abbreviated names...
        return self.record.get('journal')
        
    @property
    def volume(self):
        """
        The volume of this reference as a string.
        """
        return self.record.get('volume')
    
    @property
    def pages(self):
        """
        Get the pages of this reference as a string.
        """
        return self.record.get('pages')
    
    @property
    def sections(self):
        """
        Get a list of the sections that make up a reference of this style.
        """
        return [self.authors, self.journal, self.year, self.volume, self.pages]

class RSC_format_online(RSC_format):
    """
    Class for printing online references in RSC (Royal Society of Chemistry) format.
    """
    
    @property
    def urldate(self):
        """
        The date the resource was accessed.
        """
        if self.record.get('urldate') is None:
            return None
        
        # Parse our urldate field into an object (we probably should check the format of urldate rather than assuming it's always the same...)
        date_obj = datetime.strptime(self.record.get('urldate'), '%Y-%m-%d')
        # And return just the month and year, which is all RSC asks for.
        return date_obj.strftime("%B %Y")
        #return "{} {}".format(date_obj.month, date_obj.year)
        #return self.record.get('urldate')
    
    @property
    def sections(self):
        """
        Get a list of the sections that make up a reference of this style.
        """
        return [self.authors, self.url, self.year, "(accessed {})".format(self.urldate) if self.urldate is not None else None]

class RSC_format_thesis(RSC_format):
    """
    Class for printing thesis (PhD, Masters etc) references in RSC (Royal Society of Chemistry) format.
    """
    
    @property
    def thesis_type(self):
        """
        The type of this thesis as a string.
        
        Examples are "PhD", "Masters"
        """
        if self.record.get('ENTRYTYPE') == "mastersthesis":
            return "Masters Thesis"
        elif self.record.get('ENTRYTYPE') == "phdthesis":
            return "PhD Thesis"
        else:
            return "Thesis"
        
    @property
    def school(self):
        """
        The school of this thesis.
        """
        return self.record.get("school")
        
        
    @property
    def sections(self):
        """
        Get a list of the sections that make up a reference of this style.
        """
        return [self.authors, self.thesis_type, self.school, self.year]
        
        
class RSC_format_software(RSC_format):
    """
    Classs for printing software references in RSC (Royal Society of Chemistry) format.
    """
    
    @property
    def title(self):
        """
        The title of this software (its name).
        """
        return self.record.get("title")
    
    @property
    def publisher(self):
        """
        The publisher of this software.
        """
        return self.record.get("publisher")
    
    
    @property
    def sections(self):
        """
        Get a list of the sections that make up a reference of this style.
        """
        return [self.authors, self.title, self.publisher, self.year]
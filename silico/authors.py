# Project authorship information.

authorship = {
    "Lead": "Oliver S. Lee",
    "Supervision": "Eli Zysman-Colman",
    "Beta Testers": [
        "Campbell F. R. Mackenzie",
        "Tomas Matulaitis",
        "Ettore Crovini",
        ],
    # Additional authors to go here hopefully one day...
    "Contributors": []
    }

def get_authorship_string():
    """
    """
    auth_string = """Lead: {}
Beta Testers:
    {}
Supervision: {}
""".format(authorship["Lead"], "\n    ".join(authorship["Beta Testers"]), authorship["Supervision"])
    
    if len(authorship['Contributors']) > 0:
        auth_string += "Contributors:\n    " + "\n    ".join(authorship['Contributors'])
    
    return auth_string
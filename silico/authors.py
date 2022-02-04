# Project authorship information.

authorship = {
    "Lead": "Oliver S. Lee",
    "Supervision": "Eli Zysman-Colman",
    "Beta Testers": [
        "Ettore Crovini",
        "Campbell Mackenzie",
        "Tomas Matulaitis"
        ],
    # Additional authors to go here hopefully one day...
    "Contributors": []
    }

def get_authorship_string():
    """
    """
    auth_string = """Lead: {}
Supervision: {}
Beta Testers:
    {}
""".format(authorship["Lead"], authorship["Supervision"], "\n    ".join(authorship["Beta Testers"]))
    
    if len(authorship['Contributors']) > 0:
        auth_string += "Contributors:\n    " + "\n    ".join(authorship['Contributors'])
    
    return auth_string
## -*- coding: utf-8 -*-

## Template for displaying journal references in RSC style.

<%page args="reference, reference_title = None"/>

%if reference.authors is not None:
<div class="reference__authors">${reference.authors},</div>
%endif
%if reference.journal is not None:
<div class="reference__journal">${reference.journal},</div>
%endif
%if reference.year is not None:
<div class="reference__year">${reference.year},</div>
%endif
%if reference.volume is not None:
<div class="reference__volume">${reference.volume},</div>
%endif
%if reference.pages is not None:
<div class="reference__pages">${reference.pages}</div>
%endif

## -*- coding: utf-8 -*-

<%page args="result" />

%if not result.metadata.success or result.metadata.optimisation_converged == False:
<div class="section section--fullPage">
    <div class="section__body section__body--generalWarning">
        <div class="title title--generalWarning">
        	<h1 class="title__superTitle title__superTitle--report title__superTitle--bad">Warning: Preliminary Results Only</h1>
        	%if not result.metadata.success:
        	<h2 class="title__mainTitle title__mainTitle--report">This calculation did not complete successfully!</h2>
        	%endif
        	%if result.metadata.optimisation_converged == False:
        	<h2 class="title__mainTitle title__mainTitle--report">This optimisation did not converge!</h2>
        	%endif
        </div>
        ##<div class="generalWarning">
        ##	The results contained within may be inaccurate, nonsensical and/or misleading.
        ##</div>
    </div>
</div>
%endif
function plot_coverages(sol)
%plots fractional coverages

figure
%plot theta vs. t if QSS was not assumed
if isfield(sol,'t') == true
    plot(sol.t,sol.theta)
    ylim([0 1])
    xlabel('Time (s)')
    ylabel('Fractional coverage')
    legend(sol.site_species{:},'Location','Best')
    legend Boxoff
    
    %plot bar graph of coverages if QSS was assumed
else
    bar(sol.theta)
    ylim([0 1])
    set(gca,'xticklabel',sol.site_species);
    ylabel('Fractional coverage')
end

end
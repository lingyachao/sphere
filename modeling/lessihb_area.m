function lessihb_idx = lessihb_area

    load('./computed_sphere_grid/N10242_R10.mat', 'locs');
    
    boundary_idx = [277 299 321 325 348 349 366 379 406 425 430 447 467 470 482 513 536 568 589 609 632 685 701 710 765 767 771 778 805 829 853 862 869 871 911 930 957 963 973 984 991 1004 1009 1027 1043 1048 1066 1076 1079 1084 1133 1165 1177 1184 1208 1256 1273 1286 1309 1315 1361 1384 1398 1456 1499 1513 1526 1531 1577 1640 1641 1654 1660 1723 1741 1772 1809 1812 1880 1900 1901 1949 1962 1964 2020 2021 2037 2040 2063 2067 2116 2121 2132 2133 2182 2192 2194 2204 2216 2225 2229 2236 2237 2264 2284 2308 2316 2329 2336 2337 2349 2372 2388 2404 2412 2433 2434 2437 2484 2492 2493 2512 2540 2569 2574 2622 2623 2638 2674 2675 2679 2707 2721 2725 2759 2783 2791 2823 2851];
    k = boundary(locs(boundary_idx,1), locs(boundary_idx,2), 0.8);
    boundary_idx = boundary_idx(k);
    
    p = locs(:,3) < 0;
    
    scatter(locs(p,1), locs(p,2), 15, 'y', 'filled');
    hold on;
    plot(locs(boundary_idx,1), locs(boundary_idx,2)');
    
    [in, on] = inpolygon(locs(:,1), locs(:,2), locs(boundary_idx,1), locs(boundary_idx,2));
    lessihb_idx = (in | on) & p;
    scatter(locs(lessihb_idx,1), locs(lessihb_idx,2), 15, 'g', 'filled');
    
end


isinstalled(pkg::AbstractString) = pkg != "METADATA" && pkg != "REQUIRE" && pkg[1] != '.' && Pkg.cd(isdir, pkg)

print("\nHDF5\n")
if isinstalled("HDF5")
    print("HDF5 already installed\n")
else
    print("\nAdding HDF5 package to julia\n")
    Pkg.add("HDF5")
    Pkg.update("HDF5")
end


print("\ProgressMeter\n")
if isinstalled("ProgressMeter")
    print("ProgressMeter already installed\n")
else
    print("\nAdding ProgressMeter package to julia\n")
    Pkg.add("ProgressMeter")
    Pkg.update("ProgressMeter")
end



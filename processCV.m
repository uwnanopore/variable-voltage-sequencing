function generated_read = processCV(read)
% runs each part of the DC sequencing pipeline on a read, starting only
% from measured ion current and open state current data.

read = findLevelsCV(read);
read = removalFilterCV(read);
read = recombinationFilterCV(read);
generated_read = sequenceCV(read);

end
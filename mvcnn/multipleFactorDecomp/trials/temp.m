load('AVLetter_Train_SVD-merge_letters_reps.mat');
figure
bar(u{1})
title('Letter Vector')
figure
bar(u{2})
title('Person Style Vector')

V1 = u{1}(:,1:3);
%plotting(V1,trainLetters);
V1 = u{2}(:,1:3);
plotting(V1,trainSubjects);

return;
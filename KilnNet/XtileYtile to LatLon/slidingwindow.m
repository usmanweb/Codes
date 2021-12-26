function [] = slidingwindowmat(zoom_17_image,a,b,c)
d=c;
i=1;
j=1;
count = 1;
mg = size(zoom_17_image,1);
sb=b;

 while(j+31 <= size(zoom_17_image,1))
     i = 1;
     b =sb;
     while (i+31 <= size(zoom_17_image,2))
         if (1)   
         imwrite(imresize(zoom_17_image(j:j+31,i:i+31,:),8),fullfile('Zoom20Conversion/Nepal/WholeKathmandu/',strcat(num2str(b,16),'_',num2str(a,16),'.jpg')));
         end
         b = b + c;
         i = i + 32;
         count = count + 1;
     end
     a = a -d;
     j = j + 32;
 end
end 
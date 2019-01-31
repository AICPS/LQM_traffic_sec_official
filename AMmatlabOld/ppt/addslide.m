function addslide(op,titletext,pic_H, pic_W,neworlast,pos,prnopt)
% addslide(op,titletext,neworlast,pos,prnopt)

% Capture current figure/model into clipboard:
if nargin<7
    print -dmeta
else
    print('-dmeta',prnopt)
end

if nargin<6
    pos = [0.5 0.5];
end

if nargin<5
    neworlast = 'new';
end

% Get current number of slides:
slide_count = get(op.get('Slides'),'Count');

% Add a new slide (with title object):
if strcmp(neworlast,'new')
    slide_count = int32(double(slide_count)+1);
    slide = invoke(op.get('Slides'),'Add',slide_count,11);
else
    slide_count = int32(double(slide_count));
    slide=invoke(invoke(op.get('Slides'),'Range'),'Item',slide_count);
end

% Insert text into the title object:
if(~isempty(titletext) && strcmp(neworlast,'new') )
    set(slide.Shapes.Title.TextFrame.TextRange,'Text',titletext);
end

% Get height and width of slide:
slide_H = op.get('PageSetup').get('SlideHeight');
slide_W = op.get('PageSetup').get('SlideWidth');

% Paste the contents of the Clipboard:
pic1 = invoke(slide.get('Shapes'),'Paste');

set(pic1,'Height',pic_H)
set(pic1,'Width',pic_W)

% % Center picture on page (below title area):
set(pic1,'Left',single(max([pos(1)*double(slide_W)-0.5*double(pic_W) 0])));
set(pic1,'Top' ,single(max([(1-pos(2))*double(slide_H)-0.65*double(pic_H) 0])));


return
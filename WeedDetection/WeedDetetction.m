function varargout = WeedDetetction(varargin)
% WEEDDETETCTION MATLAB code for WeedDetetction.fig
%      WEEDDETETCTION, by itself, creates a new WEEDDETETCTION or raises the existing
%      singleton*.
%
%      H = WEEDDETETCTION returns the handle to a new WEEDDETETCTION or the handle to
%      the existing singleton*.
%
%      WEEDDETETCTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WEEDDETETCTION.M with the given input arguments.
%
%      WEEDDETETCTION('Property','Value',...) creates a new WEEDDETETCTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before WeedDetetction_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to WeedDetetction_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WeedDetetction

% Last Modified by GUIDE v2.5 23-Apr-2016 21:04:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @WeedDetetction_OpeningFcn, ...
                   'gui_OutputFcn',  @WeedDetetction_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before WeedDetetction is made visible.
function WeedDetetction_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WeedDetetction (see VARARGIN)

% Choose default command line output for WeedDetetction
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes WeedDetetction wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = WeedDetetction_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ReadImage.
function ReadImage_Callback(hObject, eventdata, handles)
% hObject    handle to ReadImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global im

[file path]=uigetfile({'*.jpg';'*.bmp';'*.jpeg';'*.png';'.tif'}, 'Load Image File within Avilable Extensions');
image=[path file];
handles.file=image;
if (file==0)
    warndlg('You did not selected any file ') ; % fille is not selected
end
[fpath,fname,fext]=fileparts(file);
 validex=({'.bmp','.jpg','.jpeg','.png','.tif'});
 found=0;
 for (x=1:length(validex))
 if (strcmpi(fext,validex{x}))
     found=1;
im=imread(image);
axes(handles.Image);cla;
% axes(handles.axes4);

h = waitbar(0,'Please wait.Image is Uploadnig...');
steps = 50;
for step = 1:steps
    % computations take place here
    waitbar(step / steps)
end
close(h) 
  end
 end
 
 imshow(im);
 axes(handles.Image);
if (found==0)
     errordlg('Selected file does not match available extensions. Please select file from available extensions [ .jpg, .jpeg, .bmp, .png] ','Image Format Error');
end

% GrayImage_Callback(hObject, eventdata, handles)
% Enhancement_Callback(hObject, eventdata, handles)
% 
% Binarization_Callback(hObject, eventdata, handles)
% AreaThresholding_Callback(hObject, eventdata, handles)

% --- Executes on button press in GrayImage.
function GrayImage_Callback(hObject, eventdata, handles)
% hObject    handle to GrayImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global im
global Igray

fontSize=15;

[rows, columns ,numberOfColorChannels] = size(im)
if numberOfColorChannels > 1
  Igray = rgb2gray(im);
else
  Igray = (im); % It's already gray.
end

axes(handles.Image);
imshow(Igray);

% h=msgbox('Operation Complete');
title('Gray Scale Image', 'FontSize', fontSize);

% --- Executes on button press in Binarization.
function Binarization_Callback(hObject, eventdata, handles)
% hObject    handle to Binarization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
global BW
inimage = im; %// type uint8
mask = inimage(:,:,1)<inimage(:,:,2) & inimage(:,:,3)<inimage(:,:,2); %// 2D mask
outimage = bsxfun(@times, inimage, uint8(mask)); %// apply mask replicated along 3rd dim
r=rgb2gray(outimage);
 level = graythresh(r);
BW = im2bw(r,level);
imshow(outimage);
imshow(BW);


% --- Executes on button press in Labeling.
function Labeling_Callback(hObject, eventdata, handles)
% hObject    handle to Labeling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Filtering.
function Filtering_Callback(hObject, eventdata, handles)
% hObject    handle to Filtering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global BW
global BW2
 K = medfilt2(BW);
BW2 = bwareaopen(K, 350);
imshow(BW2);

% --- Executes on button press in AreaThresholding.
function AreaThresholding_Callback(hObject, eventdata, handles)
% hObject    handle to AreaThresholding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global BW
global Iout

LB = 10;
% UB = 3539800099123456123;
UB = 1245453;

Iout = xor(bwareaopen(BW,LB),  bwareaopen(BW,UB));
imshow(Iout);


% imshow(J)
% BB=bwareafilt(Iout,[4 5]);
% % % % % mask=false(size(Iout));
% % % % % mask(25:end-25,25:end-25)=true;
% % % % % figure,imshow(mask);
% % % % % 
% % % % % bw=activecontour(Iout,mask,300);
% % % % % %  visboundaries(mask,'Color','r');
% % % % %  figure,imshow(bw)
%  bw3=imclearborder(BB);
% imshow(bw3);
% [r c]=size(BB);
% % z=handles.file;
% nn=imread('red.jpg');
% z=imfuse(BB,nn);
% % z=imresize(BB,[r c]);
% % nn=imresize(nn,[r c]);
% % 
% % for i=1:r
% %     for j=1:c
% %         if BB(i,j)>=254;
% %             z(i,j)=nn(i,j);
% %         else
% %             z(i,j)=z(i,j);
% %         end
% %     end
% % end
% imshow(z);
% LB = 5000;
% UB = 9000;
% IL = bwlabel(BW);
% R = regionprops(BW,'Area');
% ind = find([R.Area] >= LB & [R.Area] <= UB);
% Iout = ismember(IL,ind);
% imshow(Iout);


% BW2 = bwareaopen(BW,119000);
% 
% imshow(BW2);




% h=msgbox('Operation Complete');
% --- Executes on button press in Outpoutimage.
function Outpoutimage_Callback(hObject, eventdata, handles)
% hObject    handle to Outpoutimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in preprocessing.
function preprocessing_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Igray
a = double(Igray);
% h_dh = [-1 0 1;-2 0 2;-1 0 1];
% h_dv = [1 2 1;0 0 0;-1 -2 -1];
% h_hp = [-1/16 -1/8 -1/16;-1/8 3/4 -1/8;-1/16 -1/8 -1/16];
% for x = 1:5
%     for y = 1:5
%         h_g(x,y) = exp(-(x^2 + y^2)/(2));
%     end
% end
% 
% DH = conv2(h_dh,a);
% DV = conv2(h_dv,a);
% HP = conv2(h_hp,a);
% G = conv2(h_g,a);
% 
% 
% axes(handles.axes2);
% imshow(DH,[]);
% 
% axes(handles.axes3);
% imshow(DV,[]);
% axes(handles.axes4);
% imshow(HP,[]);
% 
% axes(handles.axes5);
% imshow(G,[]);
% 
% handles.DH = DH;
% handles.DV = DV;
% handles.HP = HP;
% handles.G = G;
% 
% guidata(hObject, handles);


% --- Executes on button press in Cong.
function Cong_Callback(hObject, eventdata, handles)
% hObject    handle to Cong (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
G= handles.G;
 h = waitbar(0.01,'Extracting phase congruency ');
 pc = phasecong(G);
 for i=1:1000
    a=(i/1000);

 waitbar(a);
    
    
end
  close(h); 

axes(handles.axes6)
imshow(pc);
handles.pc = pc;
guidata(hObject, handles);

% --- Executes on button press in Enhancement.
function Enhancement_Callback(hObject, eventdata, handles)
% hObject    handle to Enhancement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Igray
global img
img=histeq(Igray);
% axes(handles.axes10)
 imshow(img);
%  handles.img=img;
%   guidata(hObject, handles);

% --- Executes on button press in Thresholding.
function Thresholding_Callback(hObject, eventdata, handles)
% hObject    handle to Thresholding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





% --- Executes on button press in Binirayimage.
function Binirayimage_Callback(hObject, eventdata, handles)
% hObject    handle to Binirayimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
global BW
global img
inimage = im; %// type uint8
mask = inimage(:,:,1)<inimage(:,:,2) & inimage(:,:,3)<inimage(:,:,2); %// 2D mask
outimage = bsxfun(@times, inimage, uint8(mask)); %// apply mask replicated along 3rd dim
r=rgb2gray(outimage);
 level = graythresh(r);
BW = im2bw(r,level);
imshow(outimage);
imshow(BW);


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Iout
global im
J=Iout;
[B,L] = bwboundaries(J, 'noholes');
figure; imshow(J); hold on;
for k = 1:length(B),
   boundary = B{k};
   plot(boundary(:,2),boundary(:,1),'r','LineWidth',2);
end

Totalareaofweed=bwarea(Iout)
[rows, columns ,numberOfColorChannels] = size(im);
Totalareaofimage=rows*columns

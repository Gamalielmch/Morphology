function varargout = RadonShape(varargin)
% RADONSHAPE MATLAB code for RadonShape.fig
%      RADONSHAPE, by itself, creates a new RADONSHAPE or raises the existing
%      singleton*.
%
%      H = RADONSHAPE returns the handle to a new RADONSHAPE or the handle to
%      the existing singleton*.
%
%      RADONSHAPE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RADONSHAPE.M with the given input arguments.
%
%      RADONSHAPE('Property','Value',...) creates a new RADONSHAPE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RadonShape_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RadonShape_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RadonShape

% Last Modified by GUIDE v2.5 22-Oct-2019 12:18:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RadonShape_OpeningFcn, ...
                   'gui_OutputFcn',  @RadonShape_OutputFcn, ...
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


% --- Executes just before RadonShape is made visible.
function RadonShape_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RadonShape (see VARARGIN)

% Choose default command line output for RadonShape
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RadonShape wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global  menu wc um fw1 fw2 te1 te2 bt1 ru medias cova res assi cui or procces clas anno

clas=0;
procces=0; 
cui=1;
or=[8,6];
um=[40,12];
wc=[5,14];
la={'1 = Angular,     2 = Sub-angular,    3 = Sub-rounded,     4 = Rounded,     5 = Well-rounded'};
anno=annotation('textbox',[0.15,0,0.7,0.05],'String',la,'Color',[0 1 0]);
b = annotation('textbox', get(anno,'Position'));
set(anno, 'EdgeColor',[0 0 0] );
set(b, 'BackgroundColor', [0 0 0]);
set(b, 'FaceAlpha', 0.2);
uistack(anno,'top')
set(anno,'visible','off')

m_file = uimenu(hObject,'Label','File','ForegroundColor',[0 0 0]);
menu.open=uimenu(m_file,'Label','Open binary image',...
    'Callback',{@load_Callback,handles});

set(menu.open,	'UserData', 'iopen.png') 

menu.radon = uimenu(hObject,'Label','Radon','ForegroundColor',[0 0 0]);


menu.apply = uimenu(menu.radon,'Label','Apply',...
    'Enable','off','Callback',{@Radon_Callback,handles});


menu.clasi = uimenu(menu.radon,'Label','Classify',...
    'Enable','off','Separator','on','Callback',@class_Callback);

menu.setup = uimenu(menu.radon,'Label','Setup filters',...
    'Enable','off','Separator','on','Callback',{@Radon_setup_Callback,handles});

menu.setup2 = uimenu(menu.radon,'Label','Setup Models',...
    'Enable','off','Separator','on','Callback',@Radon_setupmodels_Callback);
menu.setup3 = uimenu(menu.radon,'Label','Plot Models',...
    'Enable','off','Callback',@Radon_plotmodels_Callback);

menu.expo = uimenu(hObject,'Label','Export','ForegroundColor',[0 0 0]);

menu.expdata = uimenu(menu.expo,'Label','Data',...
    'Enable','off','Callback',{@expdata_Callback,handles});
menu.expima = uimenu(menu.expo,'Label','Images',...
    'Enable','off','Callback',@expima_Callback);
menu.help = uimenu(hObject,'Label','Help','ForegroundColor',[0 0 0]);
menu.manual= uimenu(menu.help,'Label','Manual',...
    'Enable','on','Callback',@manual_Callback);

axes('Position',[0.05 0.05 0.9 0.88],'color',[0 0 0]);
axis off

set(handles.response,'visible','on','Enable','off')
set(handles.response2,'visible','on','Enable','off')


[a,~]=imread('beforei.jpg');
set(handles.beforeb,'CData',a,'Enable','off');
[a,~]=imread('nexti.jpg');
set(handles.nextb,'CData',a,'Enable','off');



fw1 = uicontrol('units','norm','pos',[0.2 0.95 0.055 0.035],'style','edit','string','5','backgroundcolor',[1 1 1],'visible','off','Enable','on','callback',@fw1_call);
fw2 = uicontrol('units','norm','pos',[0.3 0.95 0.05 0.035],'style','edit','string','14','backgroundcolor',[1 1 1],'visible','off','Enable','on','callback',@fw2_call);
te1 = uicontrol('units','norm','pos',[0.165 0.95 0.05 0.035],'style','text','string','w1','ForegroundColor',[0.3 1 0.1]);
set(te1,'BackgroundColor',[0 0 0],'Fontsize',11,'visible','off');

te2 = uicontrol('units','norm','pos',[0.26 0.95 0.05 0.035],'style','text','string','w2','ForegroundColor',[0.3 1 0.1]);
set(te2,'BackgroundColor',[0 0 0],'Fontsize',11,'visible','off');


ru = uicontrol('units','norm','pos',[0.47 0.95 0.1 0.035],'style','text','string','Running...','ForegroundColor',[0.3 1 0.1]);
set(ru,'BackgroundColor',[0 0 0],'Fontsize',11,'visible','off');


bt1 = uicontrol('units','norm','pos',[0.378 0.950 0.05 0.035],'style','pushbutton','string','ok','backgroundcolor',[1 1 1],'visible','off','Enable','on','callback',{@ok_call,handles});

medias(:,1)=[7.0559   29.6];
medias(:,2)=[8.0772  16.8];
medias(:,3)=[5.3542    4.8];
medias(:,4)=[2.649    1];
medias(:,5)=[ 0.4569        0];
cova(:,:,1)=[0.303   0.51; 0.51   58.8];
cova(:,:,2)=[2.0521    2.8049;  2.8049    9.700];
cova(:,:,3)=[2.552   2.4250; 2.4250   5.700];
cova(:,:,4)=[0.1043   -0.1462; -0.1462    1.5];
cova(:,:,5)=[0.4293   2e-11; 2e-11    0.31];
res=1;
assi=nan;
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFrame=get(handle(handles.figure1), 'javaframe');
jicon=javax.swing.ImageIcon('icon.jpg');
jFrame.setFigureIcon(jicon);

% --- Outputs from this function are returned to the command line.
function varargout = RadonShape_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function load_Callback(hObject, eventdata, handles)
global images menu assi anno ni cui procces clas 

[filenamed, pathd]= uigetfile({'*.png;*.tiff;*.tif;*.jpg;*.jpeg;*.bmp'},'Binary Image','MultiSelect', 'on');
if pathd==0
    warndlg('Image not loaded, try again');
    return 
end
set(anno,'visible','off')


if iscell(filenamed)
    ni=length(filenamed);
    images=cell(ni,1);
    for i=1:ni
        I=imread([pathd, filenamed{i}]);
        if size(I,3)>3
            I=I(:,:,1:3);
        end
        images{i}=im2bw(I);
        
    end
else
    I=imread([pathd, filenamed]);
    if size(I,3)>3
        I=I(:,:,1:3);
    end
    images{1}=im2bw(I);
    ni=1;
end

imshow(images{1})

set(handles.response,'Enable','off')
set(handles.response2,'Enable','off')
set(menu.apply,'Enable','on')
set(menu.setup,'Enable','on')
set(menu.setup2,'Enable','on')
set(menu.setup3,'Enable','on')
set(menu.clasi,'Enable','off')
set(menu.expdata,'Enable','off')
set(menu.expima,'Enable','off')
assi=nan;

if ni>1
    set(handles.nextb,'Enable','on');
    set(handles.beforeb,'Enable','on');
else
     set(handles.nextb,'Enable','off');
    set(handles.beforeb,'Enable','off');
end
cui=1;
procces=0;
clas=0;
function Radon_Callback(hObject, eventdata, handles)
global  res wc um ru menu nr anno ni cui images or imagesp procces clas

set(anno,'visible','off')
set(ru,'visible','on')
pause(0.01)

theta =0:1:180;
nr=cell(1,ni);
for i=1:ni
    nobj=regionprops(images{i},'Image','Centroid');
    Ip=false(size(images{i},1),size(images{i},2),3);
    vcen = cat(1, nobj.Centroid);
    nrt=zeros(length(nobj),2);
    for iobj=1:length(nobj)
        I=nobj(iobj).Image;
        I=normarea2(I,50000);
        I=padarray(I,[20 20],'both');
        [a,b]=size(I);
        a=max([a,b]);
        [R,~]=radon(I,theta);
        Iref=iradon(R,theta,a);
        Iref=im2bw(Iref,0.5);
        [M,N]=size(R);
        for j=1:2
            Lo=lpfilter('btw', M, N, wc(j), or(j));
            Rw=(fftshift(fft2(R))).*Lo;
            Rw=abs(ifft2(ifftshift(Rw), 'symmetric'));
            Ir=iradon(Rw,theta,a);
            Ir=im2bw(Ir,0.45);
            Idif=bwmorph(abs(Iref-Ir),'majority');
            s=regionprops(Idif,'PixelIdxList');
            %imshowpair(Ir,Iref)
            for si=1:length(s)
                if length(s(si).PixelIdxList)<um(j)
                    Idif(s(si).PixelIdxList)=0;
                end
            end
            if j==1
                nrt(iobj,j)=(sum(Idif(:))/sum(Iref(:)))*100;
            else
                nrt(iobj,j)=length(regionprops(Idif,'Image')) ;
            end
            Itemp=Ir;
            Itemp=regionprops(Itemp,'Image');
            Itemp=Itemp.Image;
            ncen=regionprops(Itemp,'Centroid');
            ncen= ncen(:).Centroid;
            xi=ceil(vcen(iobj,1)-ncen(1));
            yi=ceil(vcen(iobj,2)-ncen(2));
            xi(xi<=0)=1;
            yi(yi<=0)=1;
            Ip(yi:yi+size(Itemp,1)-1,xi:xi+size(Itemp,2)-1,j)=Itemp;

        end
        Itemp=Iref;
        Itemp=regionprops(Itemp,'Image');
        Itemp=Itemp.Image;
        ncen=regionprops(Itemp,'Centroid');
        ncen= ncen(:).Centroid;
        xi=ceil(vcen(iobj,1)-ncen(1));
        xi(xi<=0)=1;
        yi=ceil(vcen(iobj,2)-ncen(2));
        yi(yi<=0)=1;
        Ip(yi:yi+size(Itemp,1)-1,xi:xi+size(Itemp,2)-1,j+1)=Itemp;
    end
    nr(i)={nrt};
    imagesp{i}=Ip;
end
procces=1;
clas=0;
set(handles.response,'Enable','on', 'ForegroundColor',[0,1,0])
set(handles.response2,'Enable','on', 'ForegroundColor',[0,0.6,0])
res=1;
Ip=imagesp{cui};
imshowpair(Ip(:,:,3),Ip(:,:,1))
set(ru,'visible','off')
set(menu.expdata,'Enable','on')
set(menu.expima,'Enable','on')
set(menu.clasi,'Enable','on')


function Radon_setup_Callback(hObject, eventdata, handles)
global fw1 fw2 te1 te2 bt1
set (fw1,'visible','on')
set (fw2,'visible','on')
set (te1,'visible','on')
set (te2,'visible','on')
set (bt1,'visible','on')

function Radon_setupmodels_Callback(hObject, eventdata, handles)
global medias cova  clas tt yy

sal=models();
uiwait(sal)

clas=0;
if ~isempty(tt)
cova=tt;
end
if ~isempty(yy)
medias=yy;
end



function Radon_plotmodels_Callback(hObject, eventdata, handles)
global medias cova
x1max=max(medias(2,:)+4.*sqrt(squeeze(cova(2,2,:)))');
x1min=min(medias(2,:)-4.*sqrt(squeeze(cova(2,2,:)))');
x2=x1min:0.05:x1max;
x1max=max(medias(1,:)+4.*sqrt(squeeze(cova(1,1,:)))');
x1min=min(medias(1,:)-4.*sqrt(squeeze(cova(1,1,:)))');
x1=x1min:0.05:x1max;
[X1,X2] = meshgrid(x1,x2);
F=zeros(size(X1,1),size(X1,2),5);
figure('color',[1,1,1]), hold on
xlim([x1(1),x1(end)])
ylim([x2(1),x2(end)])
zlim([0,1])

ylabel({'Regions (units)'},'FontSize',10);
xlabel({'difference (%)'},'FontSize',10);
zlabel({'Density'},'FontSize',10);
grid(gca,'on');

for i=1:5
Ft = mvnpdf([X1(:) X2(:)],medias(:,i)',cova(:,:,i));
Ft = reshape(Ft,length(x2),length(x1));
F(:,:,i)=Ft;
surf(x1,x2,Ft/(max(max(Ft))))
end
view(gca,[-54.8 68.6])
M = view(gca); 
R = M(1:3,1:3); 
x = R*[1;0;0]; 
y = R*[0;1;0]; 
z = R*[0;0;1]; 
set(get(gca,'XLabel'),'rotation',360/(2*pi)*atan(x(2)/x(1))) 
set(get(gca,'YLabel'),'rotation',360/(2*pi)*atan(y(2)/y(1))) 
set(get(gca,'ZLabel'),'rotation',360/(2*pi)*atan(z(2)/z(1)))
shading interp
colormap jet
% lightangle(-45,30)
set(gcf,'Renderer','zbuffer')
set(findobj(gca,'type','surface'),...
    'FaceLighting','phong',...
    'AmbientStrength',.3,'DiffuseStrength',.8,...
    'SpecularStrength',.9,'SpecularExponent',25,...
    'BackFaceLighting','unlit')
light('Position',[-1 0 0],'Style','local')
light('Position',[-0 2 10],'Style','local')
grid(gca,'on');

% --- Executes on button press in response.
function response_Callback(hObject, eventdata, handles)
% hObject    handle to response (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global res imagesp cui clas
Ip=imagesp{cui};
res=1;
if clas==1
    class_Callback(@class_Callback, eventdata, handles)
else
imshowpair(Ip(:,:,3),Ip(:,:,1))
end
set(hObject, 'ForegroundColor',[0,1,0])
set(handles.response2, 'ForegroundColor',[0,0.6,0])

% --- Executes on button press in response2.
function response2_Callback(hObject, eventdata, handles)
% hObject    handle to response2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Ip res imagesp cui  clas
Ip=imagesp{cui}; 
res=2;
if clas==1
    class_Callback(@class_Callback, eventdata, handles)
else
imshowpair(Ip(:,:,3),Ip(:,:,2))
end
set(hObject, 'ForegroundColor',[0,1,0])
set(handles.response, 'ForegroundColor',[0,0.6,0])



function fw1_call(src,eventdata)
global wc
str=get(src,'String');
if isnan(str2double(str))
    set(src,'string',num2str(wc(1)));
    warndlg('Input must be numerical');
end

function fw2_call(src,eventdata)
global wc
str=get(src,'String');
if isnan(str2double(str))
    set(src,'string',num2str(wc(2)));
    warndlg('Input must be numerical');
end

function ok_call(hObject, eventdata, handles)
global wc fw1 fw2 ru menu  te1 te2 bt1
set (fw1,'Enable','off')
set (fw2,'Enable','off')
set (hObject,'Enable','off')
set(ru,'visible','on')
set(menu.clasi,'Enable','on')
pause(0.1)
wc(1)=str2double(get(fw1,'String'));
wc(2)=str2double(get(fw2,'String'));


set(handles.response,'visible','on', 'ForegroundColor',[0,1,0])
set(handles.response2,'visible','on', 'ForegroundColor',[0,0.6,0])
set(ru,'visible','off')
set(menu.expdata,'Enable','on')
set(menu.expima,'Enable','on')
set (fw1,'visible','off')
set (fw2,'visible','off')
set (te1,'visible','off')
set (te2,'visible','off')
set (bt1,'visible','off')
Radon_Callback(@Radon_Callback, eventdata, handles)



function expima_Callback(hObject, eventdata, handles)
global  assi  imagesp ni ru
[file, path] = uiputfile({'*.jpg';'*.png';'*.bmp';'*.tiff'});
if path==0
    return
    set(ru,'visible','off')
end
set(ru,'visible','on')
pause(0.01)
for j=1:ni
    Ip=imagesp{j};
    nobj=regionprops(Ip(:,:,1),'Centroid');
    if file~=0
        p=strfind(file,'.');
        p=p(end);
        Ic1=uint8(zeros(size(Ip,1),size(Ip,2)));
        Ic1(Ip(:,:,3)==1)=255;
        Ic2=Ic1;
        Ic3=Ic1;
        Idep=(Ip(:,:,1)-Ip(:,:,3));
        Idep=find(Idep==1);
        Ic1(Idep)=255;
        Ic2(Idep)=0;
        Ic3(Idep)=255;
        Idep=(Ip(:,:,3)-Ip(:,:,1));
        Idep=find(Idep==1);
        Ic1(Idep)=255;
        Ic2(Idep)=0;
        Ic3(Idep)=255;
        Ic1(Idep)=0;
        Ic2(Idep)=255;
        Ic3(Idep)=0;
        Ic1(:,:,2)=Ic2;
        Ic1(:,:,3)=Ic3;
        

cen=reshape([nobj(:).Centroid],2,length(nobj));
imgOut = insertInImage(uint8(Ic1), @()text(cen(1,:),cen(2,:), num2str([1:length(nobj)]')),...
     {'fontweight','bold','color','w','fontsize',26,...
 'linewidth',3,'margin',5,'edgecolor',[1 0 0],'backgroundcolor',[0.05 0.05 0.05]});
  
        

        file1=[file(1:p-1), '_image_', num2str(j),'_Resp1', file(p:end)];
        imwrite(imgOut,[path,file1]);
        Ic1=uint8(zeros(size(Ip,1),size(Ip,2)));
        Ic1(Ip(:,:,3)==1)=255;
        Ic2=Ic1;
        Ic3=Ic1;
        Idep=(Ip(:,:,2)-Ip(:,:,3));
        Idep=find(Idep==1);
        Ic1(Idep)=255;
        Ic2(Idep)=0;
        Ic3(Idep)=255;
        Idep=(Ip(:,:,3)-Ip(:,:,2));
        Idep=find(Idep==1);
        Ic1(Idep)=255;
        Ic2(Idep)=0;
        Ic3(Idep)=255;
        Ic1(Idep)=0;
        Ic2(Idep)=255;
        Ic3(Idep)=0;
        Ic1(:,:,2)=Ic2;
        Ic1(:,:,3)=Ic3;
        
        
        
cen=reshape([nobj(:).Centroid],2,length(nobj));
imgOut = insertInImage(uint8(Ic1), @()text(cen(1,:),cen(2,:), num2str([1:length(nobj)]')),...
     {'fontweight','bold','color','w','fontsize',26,...
 'linewidth',3,'margin',5,'edgecolor',[1 0 0],'backgroundcolor',[0.05 0.05 0.05]});

        file1=[file(1:p-1), '_image_', num2str(j),'_Resp2', file(p:end)];
        imwrite(imgOut,[path,file1]);
        

        
    end
end
set(ru,'visible','off')

function expdata_Callback(hObject, eventdata, handles)
global  nr assi ni images   cui2 ru
set(ru,'visible','on')
pause(0.01)

eti={'Angular';'Sub-angular';'Sub-rounded';'Rounded';'Well-rounded'};
dire = {'*.xls'};
[file, path] = uiputfile(dire);
if path==0
    set(ru,'visible','off')
    return
end
FileName=[path, file];
fileID = fopen(FileName,'w');
fprintf(fileID, ['Shape Analysis using RadonS ', date,'  \n \n']);
for j=1:ni
nobj=regionprops(images{j},'Centroid');
cen=reshape([nobj(:).Centroid],2,length(nobj));
imgOut = insertInImage(uint8(images{j}*255), @()text(cen(1,:),cen(2,:), num2str([1:length(nobj)]')),...
    {'fontweight','bold','color','w','fontsize',26,...
'linewidth',3,'margin',5,'edgecolor',[1 0 0],'backgroundcolor',[0.05 0.05 0.05]});
 imwrite(imgOut,[FileName(1:end-4),'_label_image_',num2str(j),'.jpg']);
cui2=j;
class2_Callback(@class2_Callback, eventdata, handles)
fprintf(fileID,['Image  ', num2str(j)]);
fprintf(fileID,'\n \n');
fprintf(fileID,['Particle','\t','difference (percent)','\t','Regions','\t', 'Class']);
fprintf(fileID,'\n');
nt=nr{j};
nt=[(1:size(nt,1))' nt];
for ji=1:size(nt,1)
    fprintf(fileID,'%d \t %f \t %f \t %s \n',nt(ji,1),nt(ji,2),nt(ji,3),eti{assi(ji)});
end
fprintf(fileID,'\n \n');
end
fclose(fileID);
set(ru,'visible','off')

function class_Callback(hObject, eventdata, handles)
global medias cova nr anno assi c imagesp cui clas res images

clas=1;
nrt=nr{cui};
a=size(nrt,1);
assi=uint8(zeros(a,1));
F=0;

for part=1:a
    F(1) = mvnpdf([nrt(part,1) nrt(part,2)],medias(:,1)',cova(:,:,1));
    F(2) = mvnpdf([nrt(part,1) nrt(part,2)],medias(:,2)',cova(:,:,2));
    F(3) = mvnpdf([nrt(part,1) nrt(part,2)],medias(:,3)',cova(:,:,3));
    F(4) = mvnpdf([nrt(part,1) nrt(part,2)],medias(:,4)',cova(:,:,4));
    F(5) = mvnpdf([nrt(part,1) nrt(part,2)],medias(:,5)',cova(:,:,5));
    [~,o]=max(F);
    assi(part)=o;
    
end
Ip=imagesp{cui};

if res==2
    imshowpair(Ip(:,:,3),Ip(:,:,2))
else
    imshowpair(Ip(:,:,3),Ip(:,:,1))
end

It=images{cui};
c=regionprops(It,'Centroid');
c = cat(1, c.Centroid);
for i=1:a
    hnd1=text(c(i,1),c(i,2),num2str(assi(i)));
    set(hnd1,'FontUnits','pixels','FontSize',12,'Color',[1 0 0])
end
set(anno,'visible','on')


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nextb.
function nextb_Callback(hObject, eventdata, handles)
% hObject    handle to nextb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global images ni cui res imagesp procces clas 
if cui+1<=ni
    cui=cui+1;
end

if procces==0
    imshow(images{cui})
elseif clas==1
    class_Callback(@class_Callback, eventdata, handles)
else
    Ip=imagesp{cui};
    if res==2
        imshowpair(Ip(:,:,3),Ip(:,:,2))
        set(hObject, 'ForegroundColor',[0,1,0])
        set(handles.response, 'ForegroundColor',[0,0.6,0])
    else
        imshowpair(Ip(:,:,3),Ip(:,:,1))
        set(hObject, 'ForegroundColor',[0,1,0])
        set(handles.response2, 'ForegroundColor',[0,0.6,0])
    end
end


% --- Executes on button press in beforeb.
function beforeb_Callback(hObject, eventdata, handles)
% hObject    handle to beforeb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global images cui res imagesp procces clas 

if cui-1>0
    cui=cui-1;
end
if procces==0
    imshow(images{cui})
elseif clas==1
    class_Callback(@class_Callback, eventdata, handles)
else
    Ip=imagesp{cui};
    if res==2
        imshowpair(Ip(:,:,3),Ip(:,:,2))
        set(hObject, 'ForegroundColor',[0,1,0])
        set(handles.response, 'ForegroundColor',[0,0.6,0])
    else
        imshowpair(Ip(:,:,3),Ip(:,:,1))
        set(hObject, 'ForegroundColor',[0,1,0])
        set(handles.response2, 'ForegroundColor',[0,0.6,0])
    end
    
end

function manual_Callback (hObject, eventdata, handles)

open('userguide.pdf')



function class2_Callback(hObject, eventdata, handles)
global medias cova nr anno assi c imagesp clas res cui2

clas=1;
nrt=nr{cui2};
a=size(nrt,1);
assi=uint8(zeros(a,1));
F=0;

for part=1:a
    F(1) = mvnpdf([nrt(part,1) nrt(part,2)],medias(:,1)',cova(:,:,1));
    F(2) = mvnpdf([nrt(part,1) nrt(part,2)],medias(:,2)',cova(:,:,2));
    F(3) = mvnpdf([nrt(part,1) nrt(part,2)],medias(:,3)',cova(:,:,3));
    F(4) = mvnpdf([nrt(part,1) nrt(part,2)],medias(:,4)',cova(:,:,4));
    F(5) = mvnpdf([nrt(part,1) nrt(part,2)],medias(:,5)',cova(:,:,5));
    [~,o]=max(F);
    assi(part)=o;
    
end

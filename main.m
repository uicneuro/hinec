function main(imgpath)
    nim = nim_read(imgpath);
    nim = nim_dt_spd(nim);
    nim = nim_eig(nim);
    nim = nim_fa(nim);
    
    save nim_processed.mat nim;
end

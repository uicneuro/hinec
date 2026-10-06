function update=nim_tracking_progress(label,total,interval)
% Client-side completed-seed counter. Parallel callers use DataQueue callbacks.
% Counts attempted seeds, including those that produce no retained streamline.
if nargin<3,interval=30;end
clock=tic;completed=0;last_time=-inf;last_count=-1;
update=@advance;advance(0,true);
 function advance(increment,force)
  completed=min(total,completed+increment);elapsed=toc(clock);
  if completed==last_count,return;end
  if ~force&&completed<total&&elapsed-last_time<interval,return;end
  rate=completed/max(elapsed,eps);
  if completed>=total,eta='0.0 min';
  elseif completed==0,eta='pending';
  else,eta=sprintf('%.1f min',(total-completed)/rate/60);end
  fprintf('[%s] %s: completed %d/%d seeds (%.1f%%), elapsed %.1f min, %.2f seeds/s, ETA %s\n', ...
   datestr(now,'HH:MM:SS'),label,completed,total,100*completed/max(total,1),elapsed/60,rate,eta);
  last_time=elapsed;last_count=completed;
 end
end

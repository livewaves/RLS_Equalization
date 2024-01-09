function s = chanel_state_sequence (ynoisy, env_factor);



y_env = env (abs (ynoisy (1, :)), env_factor);
chan_state = (y_env < ((max (y_env) + min (y_env)) / 2));
chan_state_changes = diff (chan_state);
modified_chan_state_changes = chan_state_changes;
modified_chan_state_changes (find (modified_chan_state_changes == -1)+env_factor) = -2;
modified_chan_state_changes (find (modified_chan_state_changes == -1)) = 0;
modified_chan_state_changes (find (modified_chan_state_changes == -2)) = -1;
s = abs (cumsum (modified_chan_state_changes));

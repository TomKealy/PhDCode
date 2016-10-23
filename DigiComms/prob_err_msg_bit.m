function pemb_t=prob_err_msg_bit(et,N,No_of_correctable_error_bits)
pemb_t=0; % Theoretical message bit error probability by Eq.(9.4.11)
for k=No_of_correctable_error_bits+1:N 
   pemb_t= pemb_t +k*nchoosek(N,k)*et.^k.*(1-et).^(N-k)/N;
end
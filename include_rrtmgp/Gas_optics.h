#ifndef GAS_OPTICS_H
#define GAS_OPTICS_H

#include "Optical_props.h"
#include "Array.h"

template<typename TF>
class Gas_optics : public Optical_props<TF>
{
    public:
        Gas_optics(
                std::vector<Gas_concs<TF>>& available_gases,
                Array<std::string,1>& gas_names,
                Array<int,3>& key_species,
                Array<int,2>& band2gpt,
                Array<TF,2>& band_lims_wavenum,
                Array<TF,1>& press_ref,
                TF press_ref_trop,
                Array<TF,1>& temp_ref,
                TF temp_ref_p,
                TF temp_ref_t,
                Array<TF,3>& vmr_ref,
                Array<TF,4>& kmajor,
                Array<TF,3>& kminor_lower,
                Array<TF,3>& kminor_upper,
                Array<std::string,1>& gas_minor,
                Array<std::string,1>& identifier_minor,
                Array<std::string,1>& minor_gases_lower,
                Array<std::string,1>& minor_gases_upper,
                Array<int,2>& minor_limits_gpt_lower,
                Array<int,2>& minor_limits_gpt_upper,
                Array<int,1>& minor_scales_with_density_lower,
                Array<int,1>& minor_scales_with_density_upper,
                Array<std::string,1>& scaling_gas_lower,
                Array<std::string,1>& scaling_gas_upper,
                Array<int,1>& scale_by_complement_lower,
                Array<int,1>& scale_by_complement_upper,
                Array<int,1>& kminor_start_lower,
                Array<int,1>& kminor_start_upper,
                Array<TF,2>& totplnk,
                Array<TF,4>& planck_frac,
                Array<TF,3>& rayl_lower,
                Array<TF,3>& rayl_upper) :
                    Optical_props<TF>(band_lims_wavenum, band2gpt),
                    totplnk(totplnk),
                    planck_frac(planck_frac)
        {
            // Temperature steps for Planck function interpolation.
            // Assumes that temperature minimum and max are the same for the absorption coefficient grid and the
            // Planck grid and the Planck grid is equally spaced.
            totplnk_delta = (temp_ref_max - temp_ref_min) / (totplnk.dim(1)-1);

            init_abs_coeffs(
                    available_gases,
                    gas_names, key_species,
                    band2gpt, band_lims_wavenum,
                    press_ref, temp_ref,
                    press_ref_trop, temp_ref_p, temp_ref_t,
                    vmr_ref,
                    kmajor, kminor_lower, kminor_upper,
                    gas_minor,identifier_minor,
                    minor_gases_lower, minor_gases_upper,
                    minor_limits_gpt_lower,
                    minor_limits_gpt_upper,
                    minor_scales_with_density_lower,
                    minor_scales_with_density_upper,
                    scaling_gas_lower, scaling_gas_upper,
                    scale_by_complement_lower,
                    scale_by_complement_upper,
                    kminor_start_lower,
                    kminor_start_upper,
                    rayl_lower, rayl_upper);
        }

    private:
        Array<TF,2> totplnk;
        Array<TF,4> planck_frac;
        TF totplnk_delta;
        TF temp_ref_min, temp_ref_max;
        TF pres_ref_min, pres_ref_max;

        Array<TF,1> press_ref, press_ref_log, temp_ref;

        Array<std::string,1> gas_names;

        void init_abs_coeffs(
                std::vector<Gas_concs<TF>>& available_gases,
                Array<std::string,1>& gas_names,
                Array<int,3>& key_species,
                Array<int,2>& band2gpt,
                Array<TF,2>& band_lims_wavenum,
                Array<TF,1>& press_ref,
                Array<TF,1>& temp_ref,
                TF press_ref_trop,
                TF temp_ref_p,
                TF temp_ref_t,
                Array<TF,3>& vmr_ref,
                Array<TF,4>& kmajor,
                Array<TF,3>& kminor_lower,
                Array<TF,3>& kminor_upper,
                Array<std::string,1>& gas_minor,
                Array<std::string,1>& identifier_minor,
                Array<std::string,1>& minor_gases_lower,
                Array<std::string,1>& minor_gases_upper,
                Array<int,2>& minor_limits_gpt_lower,
                Array<int,2>& minor_limits_gpt_upper,
                Array<int,1>& minor_scales_with_density_lower,
                Array<int,1>& minor_scales_with_density_upper,
                Array<std::string,1>& scaling_gas_lower,
                Array<std::string,1>& scaling_gas_upper,
                Array<int,1>& scale_by_complement_lower,
                Array<int,1>& scale_by_complement_upper,
                Array<int,1>& kminor_start_lower,
                Array<int,1>& kminor_start_upper,
                Array<TF,3>& rayl_lower,
                Array<TF,3>& rayl_upper)
        {
            // Which gases known to the gas optics are present in the host model (available_gases)?
            std::vector<std::string> gas_names_to_use;

            for (const std::string& s : gas_names.v())
            {
                auto it = std::find_if(
                        available_gases.begin(), available_gases.end(),
                        [&s](const auto& a){ return a.get_name() == s; } );
                if (it != available_gases.end())
                    gas_names_to_use.push_back(s);
            }

            // Now the number of gases is the union of those known to the k-distribution and provided
            // by the host model.
            const int n_gas = gas_names_to_use.size();
            Array<std::string,1> gas_names_this(std::move(gas_names_to_use), {n_gas});
            this->gas_names = gas_names_this;

            // Initialize the gas optics object, keeping only those gases known to the
            // gas optics and also present in the host model.
            // Add an offset to the indexing to interface the negative ranging of fortran.
            Array<TF,3> vmr_ref_red({vmr_ref.dim(1), n_gas+1, vmr_ref.dim(3)});
            vmr_ref_red.set_offsets({0, -1, 0});

            // allocate(vmr_ref_red(size(vmr_ref,dim=1),0:ngas, &
            //                      size(vmr_ref,dim=3)))

            //  Gas 0 is used in single-key species method, set to 1.0 (col_dry)
            for (int i1=1; i1<=vmr_ref_red.dim(1); ++i1)
                for (int i3=1; i3<=vmr_ref_red.dim(3); ++i3)
                    vmr_ref_red({i1,0,i3}) = vmr_ref({i1,1,i3});

            /*
            vmr_ref_red(:,0,:) = vmr_ref(:,1,:)
            do i = 1, ngas
              idx = string_loc_in_array(this%gas_names(i), gas_names)
              vmr_ref_red(:,i,:) = vmr_ref(:,idx+1,:)
            enddo
            call move_alloc(vmr_ref_red, this%vmr_ref)
            !
            ! Reduce minor arrays so variables only contain minor gases that are available
            ! Reduce size of minor Arrays
            !
            call reduce_minor_arrays(available_gases, &
                                     gas_names, &
                                     gas_minor,identifier_minor, &
                                     kminor_lower, &
                                     minor_gases_lower, &
                                     minor_limits_gpt_lower, &
                                     minor_scales_with_density_lower, &
                                     scaling_gas_lower, &
                                     scale_by_complement_lower, &
                                     kminor_start_lower, &
                                     this%kminor_lower, &
                                     minor_gases_lower_red, &
                                     this%minor_limits_gpt_lower, &
                                     this%minor_scales_with_density_lower, &
                                     scaling_gas_lower_red, &
                                     this%scale_by_complement_lower, &
                                     this%kminor_start_lower)
            call reduce_minor_arrays(available_gases, &
                                     gas_names, &
                                     gas_minor,identifier_minor,&
                                     kminor_upper, &
                                     minor_gases_upper, &
                                     minor_limits_gpt_upper, &
                                     minor_scales_with_density_upper, &
                                     scaling_gas_upper, &
                                     scale_by_complement_upper, &
                                     kminor_start_upper, &
                                     this%kminor_upper, &
                                     minor_gases_upper_red, &
                                     this%minor_limits_gpt_upper, &
                                     this%minor_scales_with_density_upper, &
                                     scaling_gas_upper_red, &
                                     this%scale_by_complement_upper, &
                                     this%kminor_start_upper)
        
            ! Arrays not reduced by the presence, or lack thereof, of a gas
            this%press_ref = press_ref
            this%temp_ref  = temp_ref
            this%kmajor    = kmajor
        
            if(allocated(rayl_lower) .neqv. allocated(rayl_upper)) then
              err_message = "rayl_lower and rayl_upper must have the same allocation status"
              return
            end if
            if (allocated(rayl_lower)) then
              allocate(this%krayl(size(rayl_lower,dim=1),size(rayl_lower,dim=2),size(rayl_lower,dim=3),2))
              this%krayl(:,:,:,1) = rayl_lower
              this%krayl(:,:,:,2) = rayl_upper
            end if
        
            ! ---- post processing ----
            ! Incoming coefficients file has units of Pa
            this%press_ref(:) = this%press_ref(:)
        
            ! creates log reference pressure
            allocate(this%press_ref_log(size(this%press_ref)))
            this%press_ref_log(:) = log(this%press_ref(:))
        
            ! log scale of reference pressure
            this%press_ref_trop_log = log(press_ref_trop)
        
            ! Get index of gas (if present) for determining col_gas
            call create_idx_minor(this%gas_names, gas_minor, identifier_minor, minor_gases_lower_red, &
              this%idx_minor_lower)
            call create_idx_minor(this%gas_names, gas_minor, identifier_minor, minor_gases_upper_red, &
              this%idx_minor_upper)
            ! Get index of gas (if present) that has special treatment in density scaling
            call create_idx_minor_scaling(this%gas_names, scaling_gas_lower_red, &
              this%idx_minor_scaling_lower)
            call create_idx_minor_scaling(this%gas_names, scaling_gas_upper_red, &
              this%idx_minor_scaling_upper)
        
            ! create flavor list
            ! Reduce (remap) key_species list; checks that all key gases are present in incoming
            call create_key_species_reduce(gas_names,this%gas_names, &
              key_species,key_species_red,key_species_present_init)
            err_message = check_key_species_present_init(gas_names,key_species_present_init)
            if(len_trim(err_message) /= 0) return
            ! create flavor list
            call create_flavor(key_species_red, this%flavor)
            ! create gpoint_flavor list
            call create_gpoint_flavor(key_species_red, this%get_gpoint_bands(), this%flavor, this%gpoint_flavor)
        
            ! minimum, maximum reference temperature, pressure -- assumes low-to-high ordering
            !   for T, high-to-low ordering for p
            this%temp_ref_min  = this%temp_ref (1)
            this%temp_ref_max  = this%temp_ref (size(this%temp_ref))
            this%press_ref_min = this%press_ref(size(this%press_ref))
            this%press_ref_max = this%press_ref(1)
        
            ! creates press_ref_log, temp_ref_delta
            this%press_ref_log_delta = (log(this%press_ref_min)-log(this%press_ref_max))/(size(this%press_ref)-1)
            this%temp_ref_delta      = (this%temp_ref_max-this%temp_ref_min)/(size(this%temp_ref)-1)
        
            ! Which species are key in one or more bands?
            !   this%flavor is an index into this%gas_names
            !
            if (allocated(this%is_key)) deallocate(this%is_key) ! Shouldn't ever happen...
            allocate(this%is_key(this%get_ngas()))
            this%is_key(:) = .False.
            do j = 1, size(this%flavor, 2)
              do i = 1, size(this%flavor, 1) ! should be 2
                if (this%flavor(i,j) /= 0) this%is_key(this%flavor(i,j)) = .true.
              end do
            end do
            */
        }
};
#endif

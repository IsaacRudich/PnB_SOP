"""
    sequencingRRAB(setting::RRBSetting, Union{ExactArc,SimpleArc}, packet)

Call the correct rrb method based on the settings

# Arguments
- `setting::RRBSetting`: Which function to call
- `arc::Union{ExactArc,SimpleArc}`: The current arc
- `packet`: Data with the necesary info
"""
function sequencingRRAB(setting::RRBSetting, arc::Union{ExactArc,SimpleArc}, packet)
    if setting == sop
        return sopRRAB(arc,packet[1],packet[2])
    end
end
